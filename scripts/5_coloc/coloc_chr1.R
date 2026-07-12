#!/usr/bin/env Rscript
# Colocalization at chr1:36,855,140-38,474,036 (hg19, LAVA block 35, GRIK3 region)
# PD (Nalls) x {MDD, anxiety, PTSD}. coloc.abf -> PP4 = P(shared causal variant).
suppressMessages({ if(!requireNamespace("coloc",quietly=TRUE)){cat("ERROR: coloc not installed\n");quit(status=2)}; library(coloc) })
W <- "/path/to/analysis/coloc_chr1"
rd <- function(f){ d<-read.table(file.path(W,f),header=TRUE,stringsAsFactors=FALSE); d[!duplicated(d$rsid),] }
pd <- rd("pd.tsv"); pd$varbeta<-pd$se^2; pd$maf<-pmin(pd$eaf,1-pd$eaf)
pd <- pd[pd$maf>0.01 & is.finite(pd$beta) & pd$varbeta>0 & is.finite(pd$varbeta),]
cat("PD region SNPs (MAF>1%):",nrow(pd),"\n")

run <- function(f,s2,N2,ztype,nm){
  tr <- rd(f)
  if(ztype){ f2<-pmin(tr$freq,1-tr$freq); Nn<-tr$neff
    tr$varbeta <- 1/(2*f2*(1-f2)*(Nn+tr$z^2)); tr$beta <- tr$z*sqrt(tr$varbeta); tr$maf<-f2
    if(is.na(N2)) N2 <- round(median(Nn,na.rm=TRUE))
  } else { tr$varbeta<-tr$se^2; tr$maf<-pmin(tr$eaf,1-tr$eaf) }
  tr <- tr[tr$maf>0.01 & is.finite(tr$beta) & tr$varbeta>0 & is.finite(tr$varbeta),]
  cm <- intersect(pd$rsid,tr$rsid)
  if(length(cm)<20){ cat(sprintf("PD x %s: only %d common SNPs - skip\n",nm,length(cm))); return(invisible()) }
  p<-pd[match(cm,pd$rsid),]; t<-tr[match(cm,tr$rsid),]
  d1<-list(beta=p$beta,varbeta=p$varbeta,snp=cm,type="cc",s=0.0697,N=482730,MAF=p$maf)
  d2<-list(beta=t$beta,varbeta=t$varbeta,snp=cm,type="cc",s=s2,N=N2,MAF=t$maf)
  r<-suppressWarnings(coloc.abf(d1,d2)); s<-r$summary
  cat(sprintf("\n=== PD x %s | nSNP=%d ===\n",nm,length(cm)))
  cat(sprintf("  PP0=%.3f PP1=%.3f PP2=%.3f  PP3(distinct)=%.3f  PP4(shared)=%.3f\n",
      s["PP.H0.abf"],s["PP.H1.abf"],s["PP.H2.abf"],s["PP.H3.abf"],s["PP.H4.abf"]))
  o<-r$results[order(-r$results$SNP.PP.H4),c("snp","SNP.PP.H4")][1:4,]
  cat("  top candidate shared SNPs:\n"); print(o,row.names=FALSE)
}
run("mdd.tsv",0.218,1639572,FALSE,"MDD")
run("anx.tsv",0.143,851685,FALSE,"anxiety")
run("ptsd.tsv",0.5,NA,TRUE,"PTSD")
cat("\n=== COLOC_CHR1_DONE ===\n")
