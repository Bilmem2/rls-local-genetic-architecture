#!/usr/bin/env Rscript
# Figure 2 - forest plot: RLS (Akcimen) vs PD (Nalls) global genetic correlation across the panel.
#
# Every plotted number is read from results/global_genetic_correlation_ldsc.csv (the LDSC output
# table shipped with this repository), so the figure cannot drift away from the reported results.
# This script fixes only the display order and the axis labels; if a trait is missing from the
# table it stops rather than silently plotting an incomplete panel.
ROOT <- "/path/to/analysis"
RG   <- file.path(ROOT, "results/global_genetic_correlation_ldsc.csv")
OUT  <- file.path(ROOT, "results/figures/forest_rg_panel.png")
dir.create(dirname(OUT), showWarnings=FALSE, recursive=TRUE)

d <- read.csv(RG, stringsAsFactors=FALSE, check.names=FALSE)

# display order (top to bottom) -> axis label; the values come from the table, not from this list
map <- rbind(
  c("insomnia",                     "Insomnia"),
  c("short sleep duration",         "Short sleep dur."),
  c("MDD",                          "Depression (MDD)"),
  c("PTSD",                         "PTSD"),
  c("anxiety",                      "Anxiety"),
  c("obstructive sleep apnoea",     "OSA"),
  c("long sleep duration",          "Long sleep dur."),
  c("daytime sleepiness",           "Daytime sleepiness"),
  c("REM sleep behaviour disorder", "RBD"),
  c("chronotype",                   "Chronotype"),
  c("schizophrenia",                "Schizophrenia"),
  c("bipolar disorder",             "Bipolar"),
  c("narcolepsy",                   "Narcolepsy"))

# The manuscript rounds half away from zero (e.g. rg 0.0055 -> 0.006 in Table 2); mirror that here
# so the annotation matches the text. The 1e-9 nudge absorbs binary representation error, which
# otherwise stores 0.0055 as 0.005499999... and rounds it down.
round_half_up <- function(x, d){ sign(x) * floor(abs(x) * 10^d + 0.5 + 1e-9) / 10^d }

pick <- function(anchor, trait, col){
  v <- d[d$anchor == anchor & d$trait == trait, col]
  if (length(v) != 1) stop(sprintf("expected exactly one '%s' x '%s' row in %s", anchor, trait, RG))
  as.numeric(v) }
grab <- function(anchor, col) vapply(map[,1], function(t) pick(anchor, t, col), numeric(1))

tr     <- map[,2]
rls    <- grab("RLS", "rg");                 rls.se <- grab("RLS", "se"); rls.p <- grab("RLS", "p")
pd     <- grab("Parkinson's disease", "rg"); pd.se  <- grab("Parkinson's disease", "se")
rls_pd <- pick("RLS", "Parkinson's disease", "rg")   # the null anchor-anchor correlation

n <- length(tr); y <- n:1; off <- 0.19
# Bonferroni over every trait pair tested against the anchor (as in the manuscript), not merely
# the subset plotted here. Counted from the table so the two cannot drift apart.
n_tests <- sum(d$anchor == "RLS")
col.rls <- "#1f6feb"; col.pd <- "#c0873a"
sig <- rls.p < 0.05/n_tests
png(OUT, width=2700, height=2050, res=300, type="cairo")
par(mar=c(5,9.5,3.2,1.5))
plot(NA, xlim=c(-0.42,0.63), ylim=c(0.5,n+0.5), yaxt="n", xlab="Genetic correlation rg (95% CI)", ylab="",
     main="Global genetic correlation (LDSC): RLS vs PD × sleep / psychiatric panel", cex.main=1.02)
abline(v=0, lty=2, col="gray55")
axis(2, at=y, labels=tr, las=1, cex.axis=0.92)
# RLS (upper offset)
segments(rls-1.96*rls.se, y+off, rls+1.96*rls.se, y+off, col=col.rls, lwd=2.2)
points(rls, y+off, pch=ifelse(sig,19,21), col=col.rls, bg="white", cex=1.15)
# PD (lower offset)
segments(pd-1.96*pd.se, y-off, pd+1.96*pd.se, y-off, col=col.pd, lwd=2.2)
points(pd, y-off, pch=21, col=col.pd, bg="white", cex=1.05)
legend("bottomright", legend=c("RLS (rls_akcimen)","PD (Nalls 2019)"), col=c(col.rls,col.pd),
       pch=c(19,21), lwd=2.2, bty="n", cex=0.92, pt.bg="white")
mtext(sprintf("Filled = Bonferroni-significant (RLS, p<0.05/%d).  RLS×PD global rg = %.3f (null) — their link is local (LAVA/TOX3), not global.",
              n_tests, round_half_up(rls_pd, 3)),
      side=1, line=3.6, cex=0.68, col="gray35")
dev.off()
cat("wrote", OUT, "\n")
