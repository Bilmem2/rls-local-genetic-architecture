#!/usr/bin/env bash
# =============================================================================
# TOX3 sensitivity: re-run LAVA at the TOX3 block with the UK Biobank LD reference instead of 1000G EUR.
#
# LAVA's documentation recommends the UK Biobank reference over 1000 Genomes for European analyses,
# because the small 1000G sample can bias local heritability and inflate type-1 error. Re-running
# the full panel would need a 15 GB download, so only the chr13-16 package is fetched here.
#
# Usage:
#   bash tox3_ukb_ld.sh
# The download is ~2 GB and the run takes a few minutes; re-running resumes an interrupted download.
# =============================================================================
set -uo pipefail

ROOT=/path/to/analysis
L="$ROOT/results/lava"
WD=/path/to/ukb_ld                       # download and extraction directory
URL="https://vu.data.surfsara.nl/index.php/s/DQqJ2Sqr49RP4xe/download"   # chr 13-16
ARCH="$WD/ukb_ld_chr13-16.download"
OUT="$L/tox3_ukb_ld_sensitivity.txt"
RS=/usr/bin/Rscript

mkdir -p "$WD"

# ---- 0. disk space ----
AVAIL=$(df -BG --output=avail "$WD" | tail -1 | tr -dc '0-9')
echo "### free disk: ${AVAIL} GB (at least 8 GB recommended) ###"
if [ "${AVAIL:-0}" -lt 8 ]; then
  echo "!! Low disk space. Remove this check to continue anyway."; exit 1
fi

# ---- 1. download (resumable) ----
if [ ! -s "$ARCH" ]; then
  echo "### downloading UKB LD chr13-16 (~2 GB) ###"
else
  echo "### resuming the existing download ###"
fi
curl -L -C - -o "$ARCH" "$URL" || { echo "!! download failed"; exit 1; }
ls -lh "$ARCH" | awk '{print "  indirilen:", $5}'

# ---- 2. ac ----
echo "### extracting the archive ###"
TYPE=$(file -b "$ARCH")
echo "  type: $TYPE"
case "$TYPE" in
  *Zip*)    unzip -o -q "$ARCH" -d "$WD" ;;
  *gzip*)   tar -xzf "$ARCH" -C "$WD" ;;
  *tar*)    tar -xf  "$ARCH" -C "$WD" ;;
  *)        echo "!! unrecognised archive type, extract manually: $ARCH"; exit 1 ;;
esac

# ---- 3. locate the chr16 files ----
BCOR=$(find "$WD" -name "*_chr16.bcor" | head -1)
INFO=$(find "$WD" -name "*_chr16.info" | head -1)
if [ -z "$BCOR" ] || [ -z "$INFO" ]; then
  echo "!! chr16 .bcor/.info not found. Extracted files:"; find "$WD" -maxdepth 3 -type f | head -20; exit 1
fi
PREFIX="${BCOR%.bcor}"
echo "### reference prefix: $PREFIX ###"
ls -lh "$BCOR" "$INFO" | awk '{print "  ", $9, $5}'

# ---- 4. run LAVA with the UKB reference ----
echo "### LAVA — TOX3 block, UKB LD ###"
"$RS" - "$PREFIX" "$L" "$OUT" <<'RSCRIPT'
suppressMessages(library(LAVA))
a <- commandArgs(trailingOnly = TRUE)
PREFIX <- a[1]; L <- a[2]; OUT <- a[3]

INFO <- file.path(L, "info.nalls_akcimen.txt")
OVL  <- file.path(L, "ovl.nalls_akcimen.txt")
LOCF <- file.path(L, "tox3.locfile")
for (f in c(INFO, OVL, LOCF)) if (!file.exists(f)) stop("missing input: ", f)

loci <- read.loci(LOCF)
cat("block:", sprintf("chr%s:%s-%s", loci$CHR[1], loci$START[1], loci$STOP[1]), "\n")

inp <- process.input(INFO, OVL, PREFIX, c("pd_nalls", "rls_akcimen"))
loc <- process.locus(loci[1, ], inp)
if (is.null(loc)) stop("block could not be processed (too few SNPs in the reference?)")
cat("SNPs used in the block:", loc$n.snps, " | components:", loc$K, "\n\n")

u <- run.univ(loc)
cat("=== univariate local heritability (UKB LD) ===\n"); print(u)
b <- run.bivar(loc)
cat("\n=== bivariate local rg (UKB LD) ===\n"); print(b)

cat("\n=== comparison with 1000G EUR ===\n")
cmp <- data.frame(
  reference = c("1000G EUR (published run)", "UK Biobank (this run)"),
  rho       = c(-0.5298, round(b$rho[1], 4)),
  rho.lower = c(-0.8259, round(b$rho.lower[1], 4)),
  rho.upper = c(-0.2918, round(b$rho.upper[1], 4)),
  p         = c(5.72e-05, signif(b$p[1], 3)),
  n.snps    = c(20899, loc$n.snps)
)
print(cmp, row.names = FALSE)

sink(OUT)
cat("TOX3 (chr16:52041337-53393882) — RLS(Akcimen EUR) x PD(Nalls), LD reference sensitivity\n\n")
print(u); cat("\n"); print(b); cat("\n"); print(cmp, row.names = FALSE)
sink()
cat("\n### wrote:", OUT, "###\n")
RSCRIPT

echo
echo "### done ###"
echo "Result: $OUT"
echo "To remove the download: rm -rf $WD"
