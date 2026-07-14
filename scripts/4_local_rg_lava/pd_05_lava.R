#!/usr/bin/env Rscript
# =============================================================================
# pd_05 — genome-wide LAVA LOCAL genetic correlation, PD × RLS. (novel part)
#   For each of ~2,495 LD blocks: gate on BOTH traits having local h2
#   (univ p < 0.05/2495), then bivariate local rg. Tests whether RLS-PD share
#   signal at specific loci despite the null GLOBAL rg (as ADHD-RLS did).
#   Pairs via LAVA_PAIRS = "tag,info,ovl,t1,t2;tag2,...". Generalizes 09c.
# =============================================================================
suppressMessages(library(LAVA))
ROOT <- "/path/to/analysis"
REF  <- Sys.getenv("LAVA_REF", unset=file.path(ROOT,"data/reference/g1000_eur"))   # ancestry-matched panel for cross-ancestry
LOCF <- Sys.getenv("LAVA_LOCF", unset=file.path(ROOT,"data/reference/lava_blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"))
OUT  <- Sys.getenv("LAVA_OUT",  unset=file.path(ROOT,"results/lava/lava_pd_rls.csv"))
loci <- read.loci(LOCF); n <- nrow(loci); N_GENOME <- 2495; ut <- 0.05/N_GENOME
cat("LAVA PD×RLS — loci:", n, " | both-trait univ gate p<", signif(ut,3), "\n"); flush.console()

pairs <- strsplit(Sys.getenv("LAVA_PAIRS"), ";")[[1]]
allres <- list()
for (ps in pairs) {
  f <- strsplit(ps, ",")[[1]]; tag<-f[1]; info<-f[2]; ovl<-f[3]; t1<-f[4]; t2<-f[5]
  cat(sprintf("\n>>> PAIR [%s]  %s x %s\n", tag, t1, t2)); flush.console()
  inp <- tryCatch(process.input(info, ovl, REF, c(t1,t2)),
                  error=function(e){cat("  input error:",conditionMessage(e),"\n"); NULL})
  if (is.null(inp)) next
  np <- 0
  for (i in 1:n) {
    if (i %% 250 == 0) message(sprintf("   [%s] %d/%d (hits %d)", tag, i, n, np))
    loc <- tryCatch(process.locus(loci[i,], inp), error=function(e) NULL); if (is.null(loc)) next
    u <- tryCatch(run.univ(loc), error=function(e) NULL); if (is.null(u)) next
    if (any(u$phen==t1 & u$p < ut) && any(u$phen==t2 & u$p < ut)) {
      b <- tryCatch(run.bivar(loc, phenos=c(t1,t2)), error=function(e) NULL)
      if (!is.null(b) && nrow(b) > 0) {
        b$pair<-tag; b$loc<-as.character(loci$LOC[i]); b$chr<-loc$chr; b$start<-loc$start; b$stop<-loc$stop
        allres[[paste(tag,i)]] <- b; np <- np + 1
        cat(sprintf("   [%s] %s chr%s  rho=%.3f  p=%.2g\n", tag, loci$LOC[i], loc$chr, b$rho, b$p)); flush.console()
      }
    }
  }
  cat(sprintf("   [%s] bivariate tests (both-trait gate passed): %d\n", tag, np))
}

all <- do.call(rbind, allres)
if (!is.null(all) && nrow(all) > 0) {
  all$p_bonf <- p.adjust(all$p, "bonferroni"); all$p_BH <- p.adjust(all$p, "BH")
  write.csv(all, OUT, row.names=FALSE)
  cc <- c("pair","loc","chr","phen1","phen2","rho","rho.lower","rho.upper","p","p_bonf","p_BH")
  cat("\nTotal bivariate tests:", nrow(all), "\n")
  cat("\n=== Bonferroni-significant (p_bonf<0.05) ===\n"); s<-all[which(all$p_bonf<0.05),cc];  if(nrow(s))  print(s)  else cat("  none\n")
  cat("\n=== BH-significant (p_BH<0.05) ===\n");          sb<-all[which(all$p_BH<0.05),cc]; if(nrow(sb)) print(sb) else cat("  none\n")
  cat("\n=== Top 20 by nominal p ===\n"); print(head(all[order(all$p), cc], 20))
} else cat("\nNo locus passed the both-trait univariate gate (no bivariate tests).\n")
cat("\n=== PD_LAVA_DONE ===\n")
