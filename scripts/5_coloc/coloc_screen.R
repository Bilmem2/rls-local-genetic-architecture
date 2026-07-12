#!/usr/bin/env Rscript
# Systematic coloc screen: coloc.abf EVERY significant LAVA hit (p<1e-3). Durable output -> results/coloc_screen.tsv
suppressMessages({ if(!requireNamespace("coloc",quietly=TRUE)){cat("ERROR: coloc missing\n");quit(status=2)}; library(coloc) })
ROOT<-"/path/to/analysis"; G<-file.path(ROOT,"data/gwas")
HITS<-"/path/to/analysis/lava_sig_hits.csv"; OUT<-file.path(ROOT,"results/coloc_screen.tsv")
# trait -> file, reader, case-fraction s, z-type, and awk emitting: rsid ea oa beta se [eaf]
cfg<-list(
 rls_akcimen=list(f=file.path(G,"rls/akcimen2026_EUR_RLS.txt"),rd="cat",s=0.0547,z=FALSE,
   awk='$2==c&&$3>=lo&&$3<=hi&&$10~/^rs/{print $10,$4,$5,$6,$7}'),
 pd_nalls=list(f=file.path(G,"pd/nalls2019_PD_eur.h.tsv.gz"),rd="zcat",s=0.0697,z=FALSE,
   awk='NR>1&&$1==c&&$2>=lo&&$2<=hi&&$9~/^rs/{print $9,$3,$4,$5,$6,$7}'),
 insomnia=list(f=file.path(G,"sleep/insomnia_watanabe/Watanabe_2022/insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz"),rd="zcat",s=0.283,z=FALSE,
   awk='NR>1&&$2==c&&$3>=lo&&$3<=hi&&$15~/^rs/&&$7>0{print $15,$4,$13,log($7),$8,$20}'),
 mdd=list(f=file.path(G,"daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz"),rd="zcat",s=0.218,z=FALSE,
   awk='NR>1&&$1==c&&$3>=lo&&$3<=hi&&$2~/^rs/&&$9>0{print $2,$4,$5,log($9),$10,$7}'),
 anx=list(f=file.path(G,"psych/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz"),rd="zcat",s=0.143,z=FALSE,
   awk='NR>1&&$1==c&&$3>=lo&&$3<=hi&&$2~/^rs/&&$9>0{print $2,$4,$5,log($9),$10,$7}'),
 ptsd=list(f=file.path(G,"psych/eur_ptsd_pcs_v4_aug3_2021.vcf.gz"),rd="zcat",s=0.5,z=TRUE,
   awk='$1==c&&$3>=lo&&$3<=hi&&$2~/^rs/{print $2,$4,$5,$8,$6,$7}')
)
getreg<-function(tr,chr,lo,hi){ cf<-cfg[[tr]]; if(is.null(cf))return(NULL)
  cmd<-sprintf("%s '%s' | awk -v c=%s -v lo=%.0f -v hi=%.0f '%s'",cf$rd,cf$f,chr,lo,hi,cf$awk)
  d<-tryCatch(read.table(pipe(cmd),header=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
  if(is.null(d)||nrow(d)<10)return(NULL); d<-d[!duplicated(d[[1]]),]
  if(cf$z){ names(d)<-c("rsid","ea","oa","z","freq","neff"); f<-pmin(d$freq,1-d$freq)
    d$vb<-1/(2*f*(1-f)*(d$neff+d$z^2)); d$b<-d$z*sqrt(d$vb); d$maf<-f
  } else if(ncol(d)==6){ names(d)<-c("rsid","ea","oa","b","se","eaf"); d$vb<-d$se^2; d$maf<-pmin(d$eaf,1-d$eaf)
  } else { names(d)<-c("rsid","ea","oa","b","se"); d$vb<-d$se^2; d$maf<-NA }
  d[is.finite(d$b)&is.finite(d$vb)&d$vb>0,] }
hits<-read.csv(HITS,stringsAsFactors=FALSE); res<-list()
for(k in 1:nrow(hits)){ h<-hits[k,]; t1<-h$phen1; t2<-h$phen2
  if(is.null(cfg[[t1]])||is.null(cfg[[t2]])){ cat("skip (trait not configured):",h$pair,t1,"x",t2,"\n"); next }
  a<-getreg(t1,h$chr,h$start,h$stop); b<-getreg(t2,h$chr,h$start,h$stop); if(is.null(a)||is.null(b))next
  cm<-intersect(a$rsid,b$rsid); if(length(cm)<20)next
  A<-a[match(cm,a$rsid),]; B<-b[match(cm,b$rsid),]
  maf<-B$maf; maf[!is.finite(maf)]<-A$maf[!is.finite(maf)]; ok<-is.finite(maf)&maf>0.01; if(sum(ok)<20)next
  A<-A[ok,];B<-B[ok,];cm<-cm[ok];maf<-maf[ok]
  d1<-list(beta=A$b,varbeta=A$vb,snp=cm,type="cc",s=cfg[[t1]]$s,MAF=maf)
  d2<-list(beta=B$b,varbeta=B$vb,snp=cm,type="cc",s=cfg[[t2]]$s,MAF=maf)
  r<-tryCatch(suppressWarnings(coloc.abf(d1,d2)),error=function(e)NULL); if(is.null(r))next
  s<-r$summary; top<-r$results$snp[which.max(r$results$SNP.PP.H4)]
  res[[length(res)+1]]<-data.frame(pair=h$pair,region=sprintf("chr%s:%.0f-%.0f",h$chr,h$start,h$stop),
    t1=t1,t2=t2,lava_rho=round(h$rho,3),lava_p=signif(h$p,2),nSNP=length(cm),
    PP3=round(s["PP.H3.abf"],3),PP4=round(s["PP.H4.abf"],3),top_snp=top,
    verdict=ifelse(s["PP.H4.abf"]>0.7,"COLOC",ifelse(s["PP.H4.abf"]>0.5,"weak","no"))) }
R<-do.call(rbind,res); R<-R[order(-R$PP4),]
write.table(R,OUT,sep="\t",row.names=FALSE,quote=FALSE)
cat(sprintf("\n=== COLOC SCREEN: %d hits -> %s ===\n",nrow(R),OUT)); print(R,row.names=FALSE)
cat("\n=== COLOC_SCREEN_DONE ===\n")
