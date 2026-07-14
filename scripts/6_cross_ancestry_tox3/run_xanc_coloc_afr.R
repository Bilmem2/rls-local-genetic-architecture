#!/usr/bin/env Rscript
# Cross-ancestry (AFR) local genetic correlation at the European colocalizing loci.
#
# Trait specification mirrors the main panel analysis (scripts/4_local_rg_lava/panel_lava.R):
# binary anchors (RLS, PD) carry case/control counts, every other trait is continuous
# (NA NA). The AFR anchor counts are the same ones used by the primary cross-ancestry
# scan (scripts/4_local_rg_lava/pd_xanc_run_lava.sh): rls_afr 2,176/153,313 and pd_rizig_afr 1,488/196,430.
#
# The TOX3 job is a positive control: with this specification it must reproduce the
# primary cross-ancestry estimate (rho = -0.734, CI -1.000 to -0.419) reported in the
# manuscript and Figure 5. Reference panel (g1000_afr), locus definitions (standard LAVA
# blocks) and the identity sample-overlap assumption are identical to that scan.
#
# The African-ancestry RLS arm is small (2,176 cases), so the sleep and MDD jobs are a
# directional power check, not a definitive replication test.
suppressMessages(library(LAVA))
ROOT <- "/path/to/analysis"
REF  <- file.path(ROOT,"data/reference/g1000_afr")
LOCF <- file.path(ROOT,"data/reference/lava_blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
M    <- file.path(ROOT,"results/munged")
Ldir <- file.path(ROOT,"results/lava"); OUT <- file.path(Ldir,"xanc_coloc_afr.csv")
loci <- read.loci(LOCF)

# binary anchors get case/control counts; every other trait -> continuous (NA NA)
cc <- list(rls_afr=c(2176,153313), pd_rizig_afr=c(1488,196430))
info_row <- function(ph){ f <- file.path(M, paste0(ph, ".sumstats.gz"))
  if (!is.null(cc[[ph]])) sprintf("%s %d %d %s", ph, cc[[ph]][1], cc[[ph]][2], f)
  else sprintf("%s NA NA %s", ph, f) }

JOBS <- list(
 list(name="RLS x short-sleep @ MEIS1", a="rls_afr", t="sleepdur_afr_short", loc=272,  eur="+0.43"),
 list(name="RLS x short-sleep @ PTPRD", a="rls_afr", t="sleepdur_afr_short", loc=1381, eur="+0.51"),
 list(name="RLS x MDD @ SEMA6D",        a="rls_afr", t="mdd_afr",            loc=2055, eur="+0.76"),
 list(name="RLS x PD @ TOX3",           a="rls_afr", t="pd_rizig_afr",       loc=2134, eur="-0.53"))
res <- list()
for (J in JOBS) {
  a <- J$a; t <- J$t
  info <- file.path(Ldir, "info.xanc.txt"); ovl <- file.path(Ldir, "ovl.xanc.txt")
  writeLines(c("phenotype cases controls filename", info_row(a), info_row(t)), info)
  write.table(matrix(c(1,0,0,1), 2, dimnames=list(c(a,t), c(a,t))), ovl, quote=FALSE)
  inp <- tryCatch(process.input(info, ovl, REF, c(a,t)), error=function(e){cat("input err:",conditionMessage(e),"\n");NULL})
  if (is.null(inp)) next
  loc <- tryCatch(process.locus(loci[loci$LOC==J$loc,], inp), error=function(e){cat("locus err:",conditionMessage(e),"\n");NULL})
  if (is.null(loc)) { cat(sprintf(">>> %-28s : locus not processable (no shared SNPs / variance)\n", J$name)); next }
  u <- tryCatch(run.univ(loc), error=function(e) NULL)
  ua <- if (!is.null(u)) u$p[u$phen==a][1] else NA
  ut <- if (!is.null(u)) u$p[u$phen==t][1] else NA
  b <- tryCatch(run.bivar(loc, phenos=c(a,t)), error=function(e){cat("bivar err:",conditionMessage(e),"\n");NULL})
  cat(sprintf("\n>>> %-28s  (EUR coloc rho %s)\n", J$name, J$eur))
  cat(sprintf("    univ local-h2 p:  RLS-AFR=%.2g   trait-AFR=%.2g\n", ua, ut))
  if (!is.null(b) && nrow(b) > 0) {
    cat(sprintf("    AFR local rg = %.3f  [%.3f, %.3f]   p=%.3g\n", b$rho, b$rho.lower, b$rho.upper, b$p))
    res[[length(res)+1]] <- data.frame(locus=J$name, eur_rho=J$eur, afr_rho=round(b$rho,3),
      ci_low=round(b$rho.lower,3), ci_high=round(b$rho.upper,3), p=signif(b$p,3),
      univ_p_rls=signif(ua,2), univ_p_trait=signif(ut,2))
  } else cat("    (no bivariate estimate: one trait lacks detectable local h2 here in AFR)\n")
}
if (length(res)) { R <- do.call(rbind, res); write.table(R, OUT, sep="\t", row.names=FALSE, quote=FALSE)
  cat("\n=== SUMMARY (AFR cross-ancestry at the coloc loci) ===\n"); print(R, row.names=FALSE)
  cat("\nPOSITIVE CONTROL: TOX3 must reproduce rho = -0.734 [-1.000, -0.419] (primary scan).\n") }
cat("\n=== XANC_COLOC_AFR_DONE ===\n")
