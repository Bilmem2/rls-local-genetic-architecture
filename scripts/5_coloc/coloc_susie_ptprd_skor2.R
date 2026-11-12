#!/usr/bin/env Rscript
# =====================================================================
# coloc-SuSiE at PTPRD (RLS x short sleep) and SKOR2 (RLS x daytime sleepiness).
#
# The remaining two colocalizing loci, so that the single-causal-variant
# assumption of coloc.abf is relaxed at all five prioritized loci.
# Identical to coloc_susie_tox3_sema6d.R except for the locus list and N:
# trait readers come from coloc_susie_meis1.R, window = lead +-250 kb,
# LD = 1000 Genomes EUR Phase 3 (N=503).
# =====================================================================
suppressMessages({
  library(coloc); library(susieR)
})

SRC <- Sys.getenv("SUSIE_MEIS1_SRC", "scripts/5_coloc/coloc_susie_meis1.R")   # run from the repository root
src <- readLines(SRC, warn = FALSE)
cut <- grep("^## =+$", src)
cut <- cut[cut > 60][1]                    # start of the MEIS1-specific block
if (is.na(cut)) stop("start of the MEIS1-specific block not found")
eval(parse(text = paste(src[1:(cut - 1)], collapse = "\n")))   # -> cf, getreg, mkD

PLINK <- Sys.getenv("PLINK", "plink")   # genetics PLINK 1.9 on PATH (conda env "ldsc")
BDIR  <- file.path(ROOT, "data/ldsc_ref/1000G_EUR_Phase3_plink")
WD    <- file.path(ROOT, "results/susie_ptprd_skor2"); dir.create(WD, showWarnings = FALSE)
RES   <- file.path(ROOT, "results")
RAD   <- 250000

NTOT <- list(rls_akcimen = 676175, sleepdur_eur_short = 445966, daytime_sleepiness = 452071)

LOCI <- list(
  list(tag = "ptprd", gene = "PTPRD", chr = 9,  lead = "rs12685098", partner = "sleepdur_eur_short", plab = "short sleep"),
  list(tag = "skor2", gene = "SKOR2", chr = 18, lead = "rs9965170",  partner = "daytime_sleepiness", plab = "daytime sleepiness")
)

leadpos <- function(chr, rsid) {
  bim <- file.path(BDIR, sprintf("1000G.EUR.QC.%d.bim", chr))
  p <- as.numeric(system(sprintf("awk '$2==\"%s\"{print $4;exit}' %s", rsid, shQuote(bim)), intern = TRUE))
  if (length(p) == 0 || !is.finite(p)) stop(sprintf("%s not in the reference panel", rsid))
  p
}

ldmat <- function(snps, tag, chr, lo, hi) {
  bfile <- file.path(BDIR, sprintf("1000G.EUR.QC.%d", chr))
  writeLines(snps, file.path(WD, paste0(tag, "_snps.txt")))
  sb <- file.path(WD, paste0(tag, "_sub"))
  system2(PLINK, c("--bfile", bfile, "--extract", file.path(WD, paste0(tag, "_snps.txt")),
                   "--chr", chr, "--from-bp", lo, "--to-bp", hi,
                   "--keep-allele-order", "--make-bed", "--out", sb),
          stdout = file.path(WD, "p.log"), stderr = file.path(WD, "p.log"))
  bimf <- paste0(sb, ".bim"); if (!file.exists(bimf)) return(NULL)
  bim <- read.table(bimf, stringsAsFactors = FALSE)
  names(bim) <- c("chr", "rsid", "cm", "pos", "A1", "A2")
  system2(PLINK, c("--bfile", sb, "--keep-allele-order", "--r", "square", "bin",
                   "--out", file.path(WD, tag)),
          stdout = file.path(WD, "p.log"), stderr = file.path(WD, "p.log"))
  binf <- file.path(WD, paste0(tag, ".ld.bin")); if (!file.exists(binf)) return(NULL)
  n <- nrow(bim)
  M <- matrix(readBin(binf, "double", n = n * n, size = 8), n, n)
  M[!is.finite(M)] <- 0; diag(M) <- 1
  dimnames(M) <- list(bim$rsid, bim$rsid)
  system2(PLINK, c("--bfile", sb, "--keep-allele-order", "--freq", "--out", file.path(WD, tag)),
          stdout = file.path(WD, "p.log"), stderr = file.path(WD, "p.log"))
  frqf <- file.path(WD, paste0(tag, ".frq")); refmaf <- rep(NA_real_, n)
  if (file.exists(frqf)) {
    fq <- read.table(frqf, header = TRUE, stringsAsFactors = FALSE)
    refmaf <- fq$MAF[match(bim$rsid, fq$SNP)]; refmaf <- pmin(refmaf, 1 - refmaf)
  }
  list(LD = M, bim = bim, refmaf = refmaf)
}

fitone <- function(d, bim) {
  d <- d[match(bim$rsid, d$rsid), ]
  ea <- toupper(d$ea); oa <- toupper(d$oa); a1 <- toupper(bim$A1); a2 <- toupper(bim$A2)
  same <- ea == a1 & oa == a2; flip <- ea == a2 & oa == a1
  keep <- (same | flip) & is.finite(d$b) & is.finite(d$vb) & d$vb > 0
  keep[is.na(keep)] <- FALSE
  b <- d$b; b[flip & !is.na(flip)] <- -b[flip & !is.na(flip)]
  list(beta = b, varbeta = d$vb, maf = d$maf, keep = keep,
       nsame = sum(same, na.rm = TRUE), nflip = sum(flip, na.rm = TRUE))
}

showcs <- function(f, lab) {
  if (is.null(f)) { cat(sprintf("  %-18s SuSiE failed\n", lab)); return(0) }
  n <- if (is.null(f$sets$cs)) 0 else length(f$sets$cs)
  cat(sprintf("  %-18s %d credible set", lab, n))
  if (n > 0) cat(sprintf(" -> %s", paste(sapply(f$sets$cs, function(s)
    sprintf("%s (%d SNP)", names(f$pip)[s][which.max(f$pip[s])], length(s))), collapse = "; ")))
  cat("\n"); n
}

for (L in LOCI) {
  cat(sprintf("\n================ %s : RLS x %s ================\n", L$gene, L$plab))
  pos <- leadpos(L$chr, L$lead)
  LO <- pos - RAD; HI <- pos + RAD
  cat(sprintf("  lead %s = chr%d:%s (hg19) | window %d kb\n", L$lead, L$chr, format(pos, big.mark = ","), 2 * RAD / 1000))

  a <- getreg("rls_akcimen", L$chr, LO, HI)
  b <- getreg(L$partner, L$chr, LO, HI)
  if (is.null(a) || is.null(b)) { cat("  no data\n"); next }
  cm <- intersect(a$rsid, b$rsid)
  cat(sprintf("  RLS %d SNP, %s %d SNP, shared %d\n", nrow(a), L$plab, nrow(b), length(cm)))
  if (length(cm) < 50) { cat("  too few SNPs, skipped\n"); next }

  M <- ldmat(cm, L$tag, L$chr, LO, HI)
  if (is.null(M)) { cat("  LD matrix could not be built\n"); next }
  cat(sprintf("  matched in the reference panel: %d SNP\n", nrow(M$bim)))

  fa <- fitone(a, M$bim); fb <- fitone(b, M$bim)
  cat(sprintf("  allele alignment  RLS: same %d / flipped %d | %s: same %d / flipped %d\n",
              fa$nsame, fa$nflip, L$plab, fb$nsame, fb$nflip))

  rm0 <- M$refmaf
  ok <- fa$keep & fb$keep & is.finite(rm0) & rm0 > 0.01; ok[is.na(ok)] <- FALSE
  cat(sprintf("  aligned, reference MAF>0.01: %d SNP\n", sum(ok)))
  if (sum(ok) < 50) { cat("  too few, skipped\n"); next }

  LDs <- M$LD[ok, ok, drop = FALSE]; snps <- M$bim$rsid[ok]; ps <- M$bim$pos[ok]; RMAF <- rm0[ok]
  mk <- function(f, t) {
    D <- list(beta = f$beta[ok], varbeta = f$varbeta[ok], snp = snps, MAF = RMAF,
              LD = LDs, N = NTOT[[t]], position = ps)
    if (cf[[t]]$ty == "quant") D$type <- "quant" else { D$type <- "cc"; D$s <- cf[[t]]$s }
    D
  }
  FA <- tryCatch(runsusie(mk(fa, "rls_akcimen"), n = NTOT[["rls_akcimen"]],
                          repeat_until_convergence = FALSE, max_iter = 300),
                 error = function(e) { cat("  SuSiE ERROR (RLS):", conditionMessage(e), "\n"); NULL })
  FB <- tryCatch(runsusie(mk(fb, L$partner), n = NTOT[[L$partner]],
                          repeat_until_convergence = FALSE, max_iter = 300),
                 error = function(e) { cat("  SuSiE ERROR (", L$plab, "):", conditionMessage(e), "\n"); NULL })

  na <- showcs(FA, "RLS"); nb <- showcs(FB, L$plab)
  if (is.null(FA) || is.null(FB) || na == 0 || nb == 0) { cat("  no set-level colocalization possible\n") }
  else {
    r <- tryCatch(coloc.susie(FA, FB), error = function(e) { cat("  coloc.susie ERROR:", conditionMessage(e), "\n"); NULL })
    if (is.null(r) || is.null(r$summary) || nrow(r$summary) == 0) cat("  no shared signal found\n")
    else {
      print(r$summary[, c("idx1", "idx2", "PP.H3.abf", "PP.H4.abf", "hit1", "hit2")], row.names = FALSE)
      write.csv(r$summary, file.path(RES, sprintf("susie_%s_rls_%s.csv", L$tag, L$partner)), row.names = FALSE)
    }
  }
  # single-signal comparison: coloc.abf on the same SNP set
  ab <- tryCatch(coloc.abf(mk(fa, "rls_akcimen"), mk(fb, L$partner)),
                 error = function(e) NULL)
  if (!is.null(ab)) cat(sprintf("  [coloc.abf on the same SNP set] PP.H3=%.3f PP.H4=%.3f\n",
                                ab$summary[["PP.H3.abf"]], ab$summary[["PP.H4.abf"]]))
}
cat("\n=== SUSIE_PTPRD_SKOR2_DONE ===\n")
