#!/usr/bin/env Rscript
# LAVA bivariate P/CI ne kadar stokastik? Ayni blok, 8 tekrar, iki referans.
suppressMessages(library(LAVA))
L    <- "/path/to/analysis/results/lava"
INFO <- file.path(L, "info.nalls_akcimen.txt")
OVL  <- file.path(L, "ovl.nalls_akcimen.txt")
LOCF <- file.path(L, "tox3.locfile")
REFS <- c("1000G EUR" = "/path/to/analysis/data/reference/g1000_eur",
          "UK Biobank" = "/path/to/ukb_ld/lava-ukb-v1.1_chr16")
loci <- read.loci(LOCF)
REP <- 8

out <- list()
for (nm in names(REFS)) {
  inp <- process.input(INFO, OVL, REFS[[nm]], c("pd_nalls", "rls_akcimen"))
  loc <- process.locus(loci[1, ], inp)
  rho <- p <- lo <- hi <- numeric(REP)
  up  <- numeric(REP)
  for (i in 1:REP) {
    set.seed(i)
    u <- run.univ(loc); b <- run.bivar(loc)
    rho[i] <- b$rho[1]; p[i] <- b$p[1]; lo[i] <- b$rho.lower[1]; hi[i] <- b$rho.upper[1]
    up[i]  <- u$p[u$phen == "pd_nalls"]
  }
  cat("\n########", nm, "(", REP, "tekrar, seed 1..8 )########\n")
  cat(sprintf("  rho        : hepsi ayni mi = %s   deger = %.6f\n",
              length(unique(round(rho, 8))) == 1, rho[1]))
  cat(sprintf("  bivariate P: min=%.3g  medyan=%.3g  maks=%.3g\n", min(p), median(p), max(p)))
  cat(sprintf("  CI alt     : min=%.4f  maks=%.4f\n", min(lo), max(lo)))
  cat(sprintf("  CI ust     : min=%.4f  maks=%.4f\n", min(hi), max(hi)))
  cat(sprintf("  univ P (PD): min=%.3g  medyan=%.3g  maks=%.3g  (kapi 2.0e-5)\n",
              min(up), median(up), max(up)))
  cat(sprintf("  seed=1 degerleri -> rho=%.4f  CI[%.4f, %.4f]  P=%.3g  univP(PD)=%.3g\n",
              rho[1], lo[1], hi[1], p[1], up[1]))
  out[[nm]] <- data.frame(reference = nm, rho = round(rho[1], 4),
                          lo = round(lo[1], 4), hi = round(hi[1], 4),
                          p_seed1 = signif(p[1], 3),
                          p_min = signif(min(p), 3), p_max = signif(max(p), 3),
                          univP_pd = signif(up[1], 3))
}
cat("\n\n======== seed=1 ile raporlanabilir degerler ========\n")
print(do.call(rbind, out), row.names = FALSE)
write.csv(do.call(rbind, out),
          "/path/to/analysis/results/lava/tox3_ld_reference_comparison.csv",
          row.names = FALSE)
