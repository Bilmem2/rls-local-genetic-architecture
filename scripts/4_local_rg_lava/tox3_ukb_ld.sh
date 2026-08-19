#!/usr/bin/env bash
# =============================================================================
# TOX3 duyarlilik: LAVA'yi 1000G EUR yerine UK Biobank LD referansiyla yeniden kos.
#
# NEDEN: LAVA'nin kendi dokumantasyonu Avrupa analizleri icin 1000G yerine UKB
# referansini "strongly recommend" ediyor; 1000G'nin kucuk ornekleminin lokal
# kalitsalik tahminlerinde yanlilik ve tip-1 hata sismesi yarattigi gosterilmis.
# Butun paneli yeniden kosmak 15 GB indirme demek; burada SADECE chr16 paketini
# indirip TOX3 blogunu yeniden kosuyoruz.
#
# KENDI TERMINALINDE CALISTIR:
#   bash tox3_ukb_ld.sh
# Indirme ~2 GB, kosu birkac dakika. Yarida kalirsa tekrar calistir, indirmeyi
# kaldigi yerden surdurur.
# =============================================================================
set -uo pipefail

ROOT=/path/to/analysis
L="$ROOT/results/lava"
WD=/path/to/ukb_ld                       # indirme + acma dizini
URL="https://vu.data.surfsara.nl/index.php/s/DQqJ2Sqr49RP4xe/download"   # chr 13-16
ARCH="$WD/ukb_ld_chr13-16.download"
OUT="$L/tox3_ukb_ld_sensitivity.txt"
RS=/usr/bin/Rscript

mkdir -p "$WD"

# ---- 0. disk alani ----
AVAIL=$(df -BG --output=avail "$WD" | tail -1 | tr -dc '0-9')
echo "### bos disk: ${AVAIL} GB (en az 8 GB onerilir) ###"
if [ "${AVAIL:-0}" -lt 8 ]; then
  echo "!! Disk az gorunuyor. Yine de devam etmek istersen bu kontrolu kaldir."; exit 1
fi

# ---- 1. indir (kaldigi yerden) ----
if [ ! -s "$ARCH" ]; then
  echo "### UKB LD chr13-16 indiriliyor (~2 GB) ###"
else
  echo "### mevcut indirme surduruluyor ###"
fi
curl -L -C - -o "$ARCH" "$URL" || { echo "!! indirme basarisiz"; exit 1; }
ls -lh "$ARCH" | awk '{print "  indirilen:", $5}'

# ---- 2. ac ----
echo "### arsiv aciliyor ###"
TYPE=$(file -b "$ARCH")
echo "  tur: $TYPE"
case "$TYPE" in
  *Zip*)    unzip -o -q "$ARCH" -d "$WD" ;;
  *gzip*)   tar -xzf "$ARCH" -C "$WD" ;;
  *tar*)    tar -xf  "$ARCH" -C "$WD" ;;
  *)        echo "!! taninmayan arsiv turu, elle ac: $ARCH"; exit 1 ;;
esac

# ---- 3. chr16 dosyalarini bul ----
BCOR=$(find "$WD" -name "*_chr16.bcor" | head -1)
INFO=$(find "$WD" -name "*_chr16.info" | head -1)
if [ -z "$BCOR" ] || [ -z "$INFO" ]; then
  echo "!! chr16 .bcor/.info bulunamadi. Acilan dosyalar:"; find "$WD" -maxdepth 3 -type f | head -20; exit 1
fi
PREFIX="${BCOR%.bcor}"
echo "### referans prefix: $PREFIX ###"
ls -lh "$BCOR" "$INFO" | awk '{print "  ", $9, $5}'

# ---- 4. LAVA'yi UKB referansiyla kos ----
echo "### LAVA — TOX3 blogu, UKB LD ###"
"$RS" - "$PREFIX" "$L" "$OUT" <<'RSCRIPT'
suppressMessages(library(LAVA))
a <- commandArgs(trailingOnly = TRUE)
PREFIX <- a[1]; L <- a[2]; OUT <- a[3]

INFO <- file.path(L, "info.nalls_akcimen.txt")
OVL  <- file.path(L, "ovl.nalls_akcimen.txt")
LOCF <- file.path(L, "tox3.locfile")
for (f in c(INFO, OVL, LOCF)) if (!file.exists(f)) stop("eksik girdi: ", f)

loci <- read.loci(LOCF)
cat("blok:", sprintf("chr%s:%s-%s", loci$CHR[1], loci$START[1], loci$STOP[1]), "\n")

inp <- process.input(INFO, OVL, PREFIX, c("pd_nalls", "rls_akcimen"))
loc <- process.locus(loci[1, ], inp)
if (is.null(loc)) stop("blok islenemedi (referansta yeterli SNP yok olabilir)")
cat("blokta kullanilan SNP:", loc$n.snps, " | bilesen:", loc$K, "\n\n")

u <- run.univ(loc)
cat("=== univariate lokal kalitsalik (UKB LD) ===\n"); print(u)
b <- run.bivar(loc)
cat("\n=== bivariate lokal rg (UKB LD) ===\n"); print(b)

cat("\n=== 1000G EUR ile karsilastirma ===\n")
cmp <- data.frame(
  reference = c("1000G EUR (yayindaki)", "UK Biobank (bu kosu)"),
  rho       = c(-0.5298, round(b$rho[1], 4)),
  rho.lower = c(-0.8259, round(b$rho.lower[1], 4)),
  rho.upper = c(-0.2918, round(b$rho.upper[1], 4)),
  p         = c(5.72e-05, signif(b$p[1], 3)),
  n.snps    = c(20899, loc$n.snps)
)
print(cmp, row.names = FALSE)

sink(OUT)
cat("TOX3 (chr16:52041337-53393882) — RLS(Akcimen EUR) x PD(Nalls), LD referans duyarliligi\n\n")
print(u); cat("\n"); print(b); cat("\n"); print(cmp, row.names = FALSE)
sink()
cat("\n### yazildi:", OUT, "###\n")
RSCRIPT

echo
echo "### BITTI ###"
echo "Sonuc: $OUT"
echo "Arsivi silmek istersen: rm -rf $WD"
