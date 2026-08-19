#!/usr/bin/env Rscript
# FULL systematic coloc screen: coloc.abf on EVERY p<1e-3 LAVA hit across ALL panel pairs.
# Regenerates sig-hits from all lava_panel*.csv, configs all 14 traits. Output -> results/coloc_screen_full.tsv
suppressMessages({ if(!requireNamespace("coloc",quietly=TRUE)){cat("ERROR: coloc missing\n");quit(status=2)}; library(coloc) })
ROOT<-"/path/to/analysis"; G<-file.path(ROOT,"data/gwas"); LV<-file.path(ROOT,"results/lava")
PTHR<-as.numeric(Sys.getenv("COLOC_P","0.001")); PAT<-Sys.getenv("COLOC_PAT","lava_panel"); TAG<-Sys.getenv("COLOC_TAG","full")
OUT<-file.path(ROOT,paste0("results/coloc_screen_",TAG,".tsv"))

## ---- 1. regenerate sig hits (p<PTHR) from matching LAVA csvs ----
csvs<-list.files(LV,pattern=paste0("^",PAT,".*\\.csv$"),full.names=TRUE)
H<-do.call(rbind,lapply(csvs,function(f){ d<-tryCatch(read.csv(f,stringsAsFactors=FALSE),error=function(e)NULL)
  need<-c("phen1","phen2","rho","p","pair","chr","start","stop"); if(is.null(d)||!all(need%in%names(d)))return(NULL)
  d<-d[is.finite(d$p)&d$p<PTHR,need]; d }))
H<-H[!duplicated(paste(H$pair,H$chr,H$start)),]
cat(sprintf("LAVA files: %d | sig hits (p<1e-3): %d | pairs: %d\n",length(csvs),nrow(H),length(unique(H$pair))))

## ---- 2. trait config: file, reader, type(cc/quant), s or N, flag, awk emitting: rsid ea oa c4 c5 eaf ----
cf<-list(
 rls_akcimen=list(f="rls/akcimen2026_EUR_RLS.txt",rd="cat",ty="cc",s=0.0547,fl="beta",
   awk='$2==c&&$3>=lo&&$3<=hi&&$10~/^rs/{print $10,$4,$5,$6,$7,"NA"}'),
 pd_nalls=list(f="pd/nalls2019_PD_eur.h.tsv.gz",rd="zcat",ty="cc",s=0.0697,fl="beta",
   awk='NR>1&&$1==c&&$2>=lo&&$2<=hi&&$9~/^rs/{print $9,$3,$4,$5,$6,$7}'),
 insomnia=list(f="sleep/insomnia_watanabe/Watanabe_2022/insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz",rd="zcat",ty="cc",s=0.283,fl="or",
   awk='NR>1&&$2==c&&$3>=lo&&$3<=hi&&$15~/^rs/&&$7>0{print $15,$4,$13,$7,$8,$20}'),
 mdd=list(f="daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz",rd="zcat",ty="cc",s=0.218,fl="or",
   awk='NR>1&&$1==c&&$3>=lo&&$3<=hi&&$2~/^rs/&&$9>0{print $2,$4,$5,$9,$10,$7}'),
 anx=list(f="psych/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz",rd="zcat",ty="cc",s=0.143,fl="or",
   awk='NR>1&&$1==c&&$3>=lo&&$3<=hi&&$2~/^rs/&&$9>0{print $2,$4,$5,$9,$10,$7}'),
 ptsd=list(f="psych/eur_ptsd_pcs_v4_aug3_2021.vcf.gz",rd="zcat",ty="cc",s=0.5,fl="z",
   awk='$1==c&&$3>=lo&&$3<=hi&&$2~/^rs/{print $2,$4,$5,$8,$6,$7}'),
 osa=list(f="sleep/osa/osa_rsid.tsv.gz",rd="zcat",ty="cc",s=0.15,fl="beta",
   awk='NR>1&&$2==c&&$3>=lo&&$3<=hi&&$1~/^rs/{print $1,$4,$5,$6,$7,$8}'),
 sleepdur_eur_short=list(f="sleep/sleepdur_gelernter/eur_short_rsid2.txt.gz",rd="zcat",ty="cc",s=0.27,fl="beta_p",
   awk='NR>1&&$3==c&&$4>=lo&&$4<=hi&&$1~/^rs/{print $1,$5,$6,$9,$10,$7}'),
 sleepdur_eur_long=list(f="sleep/sleepdur_gelernter/eur_long_rsid2.txt.gz",rd="zcat",ty="cc",s=0.27,fl="beta_p",
   awk='NR>1&&$3==c&&$4>=lo&&$4<=hi&&$1~/^rs/{print $1,$5,$6,$9,$10,$7}'),
 rbd=list(f="sleep/rbd_krohn/GCST90204200_buildGRCh37.tsv",rd="cat",ty="cc",s=0.112,fl="beta",
   awk='NR>1&&$1==c&&$2>=lo&&$2<=hi&&$3~/^rs/{print $3,$4,$5,$6,$7,$8}'),
 scz=list(f="psych/scz_eur_clean.tsv.gz",rd="zcat",ty="cc",s=0.409,fl="beta",
   awk='NR>1&&$1==c&&$3>=lo&&$3<=hi&&$2~/^rs/{print $2,$4,$5,$9,$10,$7}'),
 bip=list(f="psych/bip2024_eur_noUKB_no23andMe.gz",rd="zcat",ty="cc",s=0.072,fl="or",
   awk='NR>1&&$2==c&&$3>=lo&&$3<=hi&&$1~/^rs/&&$7>0{print $1,$4,$5,$7,$8,$18}'),
 chronotype=list(f="sleep/chronotype_jones/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz",rd="zcat",ty="quant",N=449734,fl="beta",
   awk='NR>1&&$2==c&&$3>=lo&&$3<=hi&&$1~/^rs/{print $1,$4,$5,$8,$9,$6}'),
 daytime_sleepiness=list(f="sleep/daytime_saxena/Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt",rd="cat",ty="quant",N=452071,fl="beta",
   awk='NR>1&&$2==c&&$3>=lo&&$3<=hi&&$1~/^rs/{print $1,$4,$5,$8,$9,$6}'),
 narcolepsy=list(f="sleep/narcolepsy_proxy/narco_clean2.tsv.gz",rd="zcat",ty="cc",s=0.0024,fl="or_p",
   awk='NR>1&&$1==c&&$2>=lo&&$2<=hi&&$9~/^rs/&&$5>0{print $9,$3,$4,$5,$8,$7}'),
 et=list(f="movement/et_verma/GCST90475832.h.tsv.gz",rd="zcat",ty="cc",s=0.031,fl="or_p",
   awk='NR>1&&$1==c&&$2>=lo&&$2<=hi&&$9~/^rs/&&$5>0{print $9,$3,$4,$5,$8,$7}'),
 tourette=list(f="movement/tourette_yu/TS_Oct2018.gz",rd="zcat",ty="cc",s=0.337,fl="or",
   awk='NR>1&&$2==c&&$3>=lo&&$3<=hi&&$1~/^rs/&&$7>0{print $1,$4,$5,$7,$8,"NA"}')
)

getreg<-function(tr,chr,lo,hi){ cc<-cf[[tr]]; if(is.null(cc))return(NULL)
  if(tr%in%c("rls_akcimen","pd_nalls")){ lo<-max(0,lo-3e6); hi<-hi+3e6 }  # widen hg38 anchor window for build-shift robustness (rsID intersect stays at trait block)
  cmd<-sprintf("%s %s 2>/dev/null | awk -v c=%s -v lo=%.0f -v hi=%.0f '%s'",cc$rd,shQuote(file.path(G,cc$f)),chr,lo,hi,cc$awk)
  d<-tryCatch(read.table(pipe(cmd),header=FALSE,stringsAsFactors=FALSE,fill=TRUE),error=function(e)NULL)
  if(is.null(d)||nrow(d)<10)return(NULL); d<-d[!duplicated(d[[1]]),]
  names(d)<-c("rsid","ea","oa","c4","c5","eaf")
  suppressWarnings({d$c4<-as.numeric(d$c4);d$c5<-as.numeric(d$c5);d$eaf<-as.numeric(d$eaf)})
  fl<-cc$fl; qn<-function(p) qnorm(1-pmin(pmax(p,1e-300),0.9999999)/2)
  if(fl=="z"){ fq<-pmin(d$c5,1-d$c5); neff<-d$eaf; vb<-1/(2*fq*(1-fq)*(neff+d$c4^2)); b<-d$c4*sqrt(vb); maf<-fq
  } else {
    if(fl=="beta"){b<-d$c4;se<-d$c5}
    else if(fl=="or"){b<-log(d$c4);se<-d$c5}
    else if(fl=="beta_p"){b<-d$c4;se<-abs(d$c4)/qn(d$c5)}
    else if(fl=="or_p"){b<-log(d$c4);se<-abs(b)/qn(d$c5)}
    vb<-se^2; maf<-pmin(d$eaf,1-d$eaf)
  }
  d$b<-as.numeric(b); d$vb<-as.numeric(vb); d$maf<-as.numeric(maf)
  d<-d[is.finite(d$b)&is.finite(d$vb)&d$vb>0,]; if(nrow(d)<10)return(NULL); d }

mkD<-function(d,tr){ cc<-cf[[tr]]; base<-list(beta=d$b,varbeta=d$vb,snp=d$rsid,MAF=d$maf)
  if(cc$ty=="quant"){base$type<-"quant";base$N<-cc$N} else {base$type<-"cc";base$s<-cc$s}; base }


## ---- prior + GERCEK genomik pencere duyarliligi ----
## Onceki surumde pencere, satir sirasina gore kirpiliyordu; intersect() sirasi
## konum sirasi olmadigi icin bu gecersizdi. Simdi blok koordinatlari daraltiliyor.
TARGET<-data.frame(
  t1=c("rls_akcimen","rls_akcimen","rls_akcimen","rls_akcimen","rls_akcimen","rls_akcimen","rls_akcimen","pd_nalls"),
  t2=c("insomnia","chronotype","sleepdur_eur_short","mdd","pd_nalls","daytime_sleepiness","sleepdur_eur_short","rbd"),
  chr=c(2,2,2,15,16,18,9,1),
  start=c(65938003,65938003,65938003,47693317,52041337,44300146,8262304,154685546),
  gene=c("MEIS1","MEIS1","MEIS1","SEMA6D","TOX3","SKOR2","PTPRD","GBA"),
  stringsAsFactors=FALSE)

P12<-c(1e-6,5e-6,1e-5,2e-5,5e-5)
FRAC<-c(1.00,0.75,0.50)

out<-list()
for(i in 1:nrow(TARGET)){
  tg<-TARGET[i,]
  hh<-H[H$phen1==tg$t1 & H$phen2==tg$t2 & H$chr==tg$chr & abs(H$start-tg$start)<1000,]
  if(nrow(hh)==0){ cat("BULUNAMADI:",tg$gene,"\n"); next }
  h<-hh[1,]; mid<-(h$start+h$stop)/2; half<-(h$stop-h$start)/2
  for(fr in FRAC){
    lo<-round(mid-half*fr); hi<-round(mid+half*fr)          # blogu konuma gore daralt
    a<-getreg(tg$t1,h$chr,lo,hi); b<-getreg(tg$t2,h$chr,lo,hi)
    if(is.null(a)||is.null(b)) next
    cm<-intersect(a$rsid,b$rsid); if(length(cm)<20) next
    A<-a[match(cm,a$rsid),]; B<-b[match(cm,b$rsid),]
    mf<-B$maf; mf[!is.finite(mf)]<-A$maf[!is.finite(mf)]; ok<-is.finite(mf)&mf>0.01
    if(sum(ok)<20) next
    A<-A[ok,];B<-B[ok,];A$maf<-mf[ok];B$maf<-mf[ok]
    for(p12 in P12){
      r<-tryCatch(suppressWarnings(coloc.abf(mkD(A,tg$t1),mkD(B,tg$t2),p12=p12)),error=function(e)NULL)
      if(is.null(r))next
      s<-r$summary
      out[[length(out)+1]]<-data.frame(gene=tg$gene,pair=paste(tg$t1,"x",tg$t2),
        window=sprintf("%.0f%%",fr*100), span_kb=round((hi-lo)/1000), nSNP=nrow(A), p12=p12,
        PP3=round(as.numeric(s["PP.H3.abf"]),3), PP4=round(as.numeric(s["PP.H4.abf"]),3),
        top=r$results$snp[which.max(r$results$SNP.PP.H4)], stringsAsFactors=FALSE)
    }
  }
  cat("bitti:",tg$gene,tg$t2,"\n")
}
R<-do.call(rbind,out)
write.table(R,file.path(ROOT,"results/coloc_prior_sensitivity.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
cat("\n=== PRIOR (tam blok) ===\n")
print(R[R$window=="100%",c("gene","pair","p12","nSNP","PP3","PP4","top")],row.names=FALSE)
cat("\n=== GENOMIK PENCERE (p12=1e-5) ===\n")
print(R[R$p12==1e-5,c("gene","pair","window","span_kb","nSNP","PP3","PP4","top")],row.names=FALSE)
cat("\n=== HER LOKUSTA PP4 ARALIGI ===\n")
ag<-aggregate(PP4~gene+pair,R,function(x)c(min(x),max(x)))
print(data.frame(gene=ag$gene,pair=ag$pair,PP4_min=ag$PP4[,1],PP4_max=ag$PP4[,2]),row.names=FALSE)
cat("\n=== PRIOR_SENS2_DONE ===\n")
