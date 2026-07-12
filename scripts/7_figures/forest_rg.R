#!/usr/bin/env Rscript
# Forest plot: RLS (rls_akcimen) vs PD (Nalls) global genetic correlation x panel
ROOT <- "/path/to/analysis"
dir.create(file.path(ROOT,"results/figures"), showWarnings=FALSE, recursive=TRUE)
tr    <- c("Insomnia","Short sleep dur.","Depression (MDD)","PTSD","Anxiety","OSA",
           "Long sleep dur.","Daytime sleepiness","RBD","Chronotype","Schizophrenia","Bipolar","Narcolepsy")
rls   <- c(0.48,0.35,0.33,0.31,0.28,0.25,0.23,0.15,0.15,-0.05,-0.06,-0.03,0.04)
rls.se<- c(0.034,0.033,0.032,0.032,0.030,0.028,0.055,0.032,0.128,0.030,0.027,0.030,0.115)
rls.p <- c(4e-45,2e-26,8e-25,4e-22,4e-22,2e-18,2e-5,4e-6,0.25,0.13,0.025,0.40,0.73)
pd    <- c(0.025,-0.054,-0.028,-0.025,-0.001,-0.026,0.019,-0.026,0.24,-0.034,-0.007,0.048,0.012)
pd.se <- c(0.042,0.042,0.032,0.036,0.034,0.034,0.083,0.039,0.188,0.030,0.029,0.038,0.138)
n <- length(tr); y <- n:1; off <- 0.19
col.rls <- "#1f6feb"; col.pd <- "#c0873a"
sig <- rls.p < 0.05/n
png(file.path(ROOT,"results/figures/forest_rg_panel.png"), width=2700, height=2050, res=300, type="cairo")
par(mar=c(5,9.5,3.2,1.5))
plot(NA, xlim=c(-0.15,0.63), ylim=c(0.5,n+0.5), yaxt="n", xlab="Genetic correlation rg (95% CI)", ylab="",
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
mtext("Filled = Bonferroni-significant (RLS, p<0.05/13).  RLS×PD global rg = 0.006 (null) — their link is local (LAVA/TOX3), not global.",
      side=1, line=3.6, cex=0.68, col="gray35")
dev.off()
cat("wrote results/figures/forest_rg_panel.png\n")
