#!/usr/bin/env Rscript
# Coloc RLS(akcimen) x PD(Nalls) at new LAVA candidate loci + TOX3 positive control
suppressMessages({ if(!requireNamespace("coloc",quietly=TRUE)){cat("ERROR: coloc missing\n");quit(status=2)}; library(coloc) })
W <- "/path/to/analysis/coloc_rlspd"
regions <- list(
  c("chr11:133.5Mb (LAVA rho=0.71)","rls_chr11.tsv","pd_chr11.tsv"),
  c("chr12:24Mb (LAVA rho=0.70)",   "rls_chr12.tsv","pd_chr12.tsv"),
  c("chr16:TOX3 (positive control)","rls_tox3.tsv", "pd_tox3.tsv"))
rd <- function(f,cols){ d<-read.table(file.path(W,f),header=FALSE,stringsAsFactors=FALSE); names(d)<-cols; d[!duplicated(d$rsid),] }
for(rg in regions){
  nm<-rg[1]
  rls<-rd(rg[2],c("rsid","ea","oa","beta","se"))
  pd <-rd(rg[3],c("rsid","ea","oa","beta","se","eaf"))
  rls$vb<-rls$se^2; pd$vb<-pd$se^2; pd$maf<-pmin(pd$eaf,1-pd$eaf)
  cm<-intersect(rls$rsid,pd$rsid)
  if(length(cm)<20){cat("\n===",nm,": only",length(cm),"common SNPs - skip\n");next}
  r<-rls[match(cm,rls$rsid),]; p<-pd[match(cm,pd$rsid),]
  ok<-is.finite(r$beta)&r$vb>0&is.finite(p$beta)&p$vb>0&p$maf>0.01
  r<-r[ok,];p<-p[ok,];cm<-cm[ok]
  d1<-list(beta=r$beta,varbeta=r$vb,snp=cm,type="cc",s=0.0547,N=676175,MAF=p$maf)
  d2<-list(beta=p$beta,varbeta=p$vb,snp=cm,type="cc",s=0.0697,N=482730,MAF=p$maf)
  res<-suppressWarnings(coloc.abf(d1,d2)); s<-res$summary
  cat(sprintf("\n=== RLS x PD @ %s | nSNP=%d ===\n",nm,length(cm)))
  cat(sprintf("  PP0=%.3f PP1=%.3f PP2=%.3f  PP3(distinct)=%.3f  PP4(SHARED)=%.3f\n",
      s["PP.H0.abf"],s["PP.H1.abf"],s["PP.H2.abf"],s["PP.H3.abf"],s["PP.H4.abf"]))
  o<-res$results[order(-res$results$SNP.PP.H4),c("snp","SNP.PP.H4")][1:4,]
  cat("  top candidate shared SNPs:\n"); print(o,row.names=FALSE)
}
cat("\n=== COLOC_RLSPD_DONE ===\n")
