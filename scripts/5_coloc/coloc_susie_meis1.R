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


## =====================================================================
## MEIS1 coloc-SuSiE (ikili). coloc.abf blok basina TEK nedensel varyant
## varsayar; SuSiE bu varsayimi kaldirip her traitte birden fazla bagimsiz
## sinyal (credible set) arar ve set-set kolokalizasyon yapar.
## Pencere: rs113851554 +-250 kb (standart ince-haritalama penceresi).
## LD: 1000G EUR Faz 3 (N=503) -- kucuk panel, sonuclar buna gore yorumlanmali.
## =====================================================================
suppressMessages(library(susieR))
PLINK<-Sys.getenv("PLINK","plink")   # conda env "ldsc" puts genetics PLINK 1.9 on PATH
BFILE<-file.path(ROOT,"data/ldsc_ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.2")
WD<-file.path(ROOT,"results/susie_meis1"); dir.create(WD,showWarnings=FALSE)
CHR<-2; LEAD<-66750564; RAD<-250000
LO<-LEAD-RAD; HI<-LEAD+RAD
NTOT<-list(rls_akcimen=676175, insomnia=386988, chronotype=449734, sleepdur_eur_short=445966)
LABEL<-c(rls_akcimen="RLS", insomnia="insomnia", chronotype="chronotype", sleepdur_eur_short="short sleep")
PARTNERS<-c("insomnia","chronotype","sleepdur_eur_short")
cat(sprintf("pencere: chr%d:%d-%d (%d kb)\n",CHR,LO,HI,(HI-LO)/1000))

ldmat<-function(snps,tag){
  writeLines(snps,file.path(WD,paste0(tag,"_snps.txt")))
  sb<-file.path(WD,paste0(tag,"_sub"))
  system2(PLINK,c("--bfile",BFILE,"--extract",file.path(WD,paste0(tag,"_snps.txt")),
                  "--chr",CHR,"--from-bp",LO,"--to-bp",HI,
                  "--keep-allele-order","--make-bed","--out",sb),
          stdout=file.path(WD,"p.log"),stderr=file.path(WD,"p.log"))
  bimf<-paste0(sb,".bim")
  if(!file.exists(bimf)) return(NULL)
  bim<-read.table(bimf,stringsAsFactors=FALSE)
  names(bim)<-c("chr","rsid","cm","pos","A1","A2")
  system2(PLINK,c("--bfile",sb,"--keep-allele-order","--r","square","bin","--out",file.path(WD,tag)),
          stdout=file.path(WD,"p.log"),stderr=file.path(WD,"p.log"))
  binf<-file.path(WD,paste0(tag,".ld.bin"))
  if(!file.exists(binf)) return(NULL)
  n<-nrow(bim)
  M<-matrix(readBin(binf,"double",n=n*n,size=8),n,n)
  M[!is.finite(M)]<-0; diag(M)<-1
  dimnames(M)<-list(bim$rsid,bim$rsid)
  # referans panel frekansi (RLS sumstats'inda eaf yok; LD ile ayni kaynak olsun)
  system2(PLINK,c("--bfile",sb,"--keep-allele-order","--freq","--out",file.path(WD,tag)),
          stdout=file.path(WD,"p.log"),stderr=file.path(WD,"p.log"))
  frqf<-file.path(WD,paste0(tag,".frq"))
  refmaf<-rep(NA_real_,n)
  if(file.exists(frqf)){
    fq<-read.table(frqf,header=TRUE,stringsAsFactors=FALSE)
    refmaf<-fq$MAF[match(bim$rsid,fq$SNP)]
    refmaf<-pmin(refmaf,1-refmaf)
  }
  list(LD=M,bim=bim,refmaf=refmaf)
}

fitone<-function(t,d,bim,LD){
  d<-d[match(bim$rsid,d$rsid),]
  ea<-toupper(d$ea); oa<-toupper(d$oa); a1<-toupper(bim$A1); a2<-toupper(bim$A2)
  same<-ea==a1&oa==a2; flip<-ea==a2&oa==a1; keep<-(same|flip)&is.finite(d$b)&is.finite(d$vb)&d$vb>0
  keep[is.na(keep)]<-FALSE
  b<-d$b; b[flip&!is.na(flip)]<- -b[flip&!is.na(flip)]
  list(beta=b,varbeta=d$vb,maf=d$maf,keep=keep,
       nsame=sum(same,na.rm=TRUE),nflip=sum(flip,na.rm=TRUE))
}

for(pt in PARTNERS){
  cat(sprintf("\n================ RLS x %s ================\n",LABEL[[pt]]))
  a<-getreg("rls_akcimen",CHR,LO,HI); b<-getreg(pt,CHR,LO,HI)
  if(is.null(a)||is.null(b)){cat("  veri yok\n");next}
  cm<-intersect(a$rsid,b$rsid)
  cat(sprintf("  RLS %d SNP, %s %d SNP, ortak %d\n",nrow(a),LABEL[[pt]],nrow(b),length(cm)))
  if(length(cm)<50){cat("  cok az SNP, atlandi\n");next}
  L<-ldmat(cm,pt); if(is.null(L)){cat("  LD matrisi uretilemedi\n");next}
  cat(sprintf("  referans panelde eslesen: %d SNP\n",nrow(L$bim)))
  fa<-fitone("rls_akcimen",a,L$bim,L$LD); fb<-fitone(pt,b,L$bim,L$LD)
  cat(sprintf("  alel hizalama  RLS: ayni %d / ters %d | %s: ayni %d / ters %d\n",
      fa$nsame,fa$nflip,LABEL[[pt]],fb$nsame,fb$nflip))
  rm0<-L$refmaf
  ok<-fa$keep&fb$keep&is.finite(rm0)&rm0>0.01
  ok[is.na(ok)]<-FALSE
  cat(sprintf("  hizalanan + referans MAF>0.01: %d SNP\n",sum(ok)))
  if(sum(ok)<50){cat("  cok az, atlandi\n");next}
  LDs<-L$LD[ok,ok,drop=FALSE]; snps<-L$bim$rsid[ok]; pos<-L$bim$pos[ok]; RMAF<-rm0[ok]
  mk<-function(f,t){
    D<-list(beta=f$beta[ok],varbeta=f$varbeta[ok],snp=snps,MAF=RMAF,
            LD=LDs,N=NTOT[[t]],position=pos)
    if(cf[[t]]$ty=="quant"){D$type<-"quant"} else {D$type<-"cc"; D$s<-cf[[t]]$s}
    D }
  FA<-tryCatch(runsusie(mk(fa,"rls_akcimen"),n=NTOT[["rls_akcimen"]],repeat_until_convergence=FALSE,max_iter=300),
               error=function(e){cat("  SuSiE HATA (RLS):",conditionMessage(e),"\n");NULL})
  FB<-tryCatch(runsusie(mk(fb,pt),n=NTOT[[pt]],repeat_until_convergence=FALSE,max_iter=300),
               error=function(e){cat("  SuSiE HATA (",pt,"):",conditionMessage(e),"\n");NULL})
  showcs<-function(f,lab){
    if(is.null(f)){cat(sprintf("  %-12s SuSiE basarisiz\n",lab));return(0)}
    n<-if(is.null(f$sets$cs)) 0 else length(f$sets$cs)
    cat(sprintf("  %-12s %d credible set",lab,n))
    if(n>0) cat(sprintf(" -> %s",paste(sapply(f$sets$cs,function(s)
      sprintf("%s (%d SNP)",names(f$pip)[s][which.max(f$pip[s])],length(s))),collapse="; ")))
    cat("\n"); n }
  na<-showcs(FA,"RLS"); nb<-showcs(FB,LABEL[[pt]])
  if(is.null(FA)||is.null(FB)||na==0||nb==0){cat("  set-set kolokalizasyon yapilamadi\n");next}
  r<-tryCatch(coloc.susie(FA,FB),error=function(e){cat("  coloc.susie HATA:",conditionMessage(e),"\n");NULL})
  if(is.null(r)||is.null(r$summary)||nrow(r$summary)==0){cat("  paylasilan sinyal bulunamadi\n");next}
  print(r$summary[,c("idx1","idx2","PP.H3.abf","PP.H4.abf","hit1","hit2")],row.names=FALSE)
  write.csv(r$summary,sprintf(file.path(ROOT,"results/susie_meis1_rls_%s.csv"),pt),row.names=FALSE)
}
cat("\n=== SUSIE_MEIS1_DONE ===\n")
