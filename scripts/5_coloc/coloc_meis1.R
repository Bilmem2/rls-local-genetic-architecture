#!/usr/bin/env Rscript
suppressMessages({ if(!requireNamespace("coloc",quietly=TRUE)){cat("ERROR: coloc missing\n");quit(status=2)}; library(coloc) })
W<-"/path/to/analysis/coloc_meis1"
rd<-function(f,cols){d<-read.table(file.path(W,f),header=FALSE,stringsAsFactors=FALSE);names(d)<-cols;d[!duplicated(d$rsid),]}
rls<-rd("rls.tsv",c("rsid","ea","oa","beta","se")); ins<-rd("ins.tsv",c("rsid","ea","oa","beta","se","maf"))
rls$vb<-rls$se^2; ins$vb<-ins$se^2
cm<-intersect(rls$rsid,ins$rsid)
r<-rls[match(cm,rls$rsid),]; i<-ins[match(cm,ins$rsid),]
ok<-is.finite(r$beta)&r$vb>0&is.finite(i$beta)&i$vb>0&i$maf>0.01
r<-r[ok,];i<-i[ok,];cm<-cm[ok]
d1<-list(beta=r$beta,varbeta=r$vb,snp=cm,type="cc",s=0.0547,N=676175,MAF=i$maf)
d2<-list(beta=i$beta,varbeta=i$vb,snp=cm,type="cc",s=0.283,N=386988,MAF=i$maf)
res<-suppressWarnings(coloc.abf(d1,d2)); s<-res$summary
cat(sprintf("\n=== RLS x insomnia @ MEIS1 (chr2:65.9-67.2Mb) | nSNP=%d ===\n",length(cm)))
cat(sprintf("  PP0=%.3f PP1=%.3f PP2=%.3f  PP3(distinct)=%.3f  PP4(SHARED)=%.3f\n",
    s["PP.H0.abf"],s["PP.H1.abf"],s["PP.H2.abf"],s["PP.H3.abf"],s["PP.H4.abf"]))
o<-res$results[order(-res$results$SNP.PP.H4),c("snp","SNP.PP.H4")][1:5,]
cat("  top candidate shared SNPs:\n"); print(o,row.names=FALSE)
cat("=== COLOC_MEIS1_DONE ===\n")
