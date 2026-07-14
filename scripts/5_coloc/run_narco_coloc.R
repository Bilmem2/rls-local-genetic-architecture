#!/usr/bin/env Rscript
# Coloc for the narcolepsy (Ollila) LAVA hits (p<1e-3). narco EUR file has no chr column,
# so region SNPs are taken from the anchor (chr:pos) and matched into narco by rsID.
suppressMessages({ if(!requireNamespace("coloc",quietly=TRUE)){cat("ERROR: coloc missing\n");quit(status=2)}; library(coloc) })
ROOT<-"/path/to/analysis"; G<-file.path(ROOT,"data/gwas")
LAVA<-file.path(ROOT,"results/lava/lava_narco_ollila.csv")
NARCO<-file.path(G,"sleep/narcolepsy_ollila/CHR_Whites_info0.7.noHLA.noSig.hetP1e-5.LDSC.N.txt.gz")
OUT<-file.path(ROOT,"results/coloc_narco_ollila.tsv")

# anchor extractors: emit rsid beta se over a chr:pos window
acfg<-list(
 rls_akcimen=list(f=file.path(G,"rls/akcimen2026_EUR_RLS.txt"),rd="cat",s=0.0547,N=676175,
   awk='$2==c&&$3>=lo&&$3<=hi&&$10~/^rs/{print $10,$6,$7}'),
 pd_nalls=list(f=file.path(G,"pd/nalls2019_PD_eur.h.tsv.gz"),rd="zcat",s=0.0697,N=482730,
   awk='NR>1&&$1==c&&$2>=lo&&$2<=hi&&$9~/^rs/{print $9,$5,$6}'))
NARCO_S<-0.05; NARCO_N<-61727   # Ollila EUR meta: Neff 12,219 -> ~3,050 cases / ~58,700 controls (s~0.05)

get_anchor<-function(tr,chr,lo,hi){cf<-acfg[[tr]]
  cmd<-sprintf("%s '%s' | awk -v c=%s -v lo=%.0f -v hi=%.0f '%s'",cf$rd,cf$f,chr,lo,hi,cf$awk)
  d<-tryCatch(read.table(pipe(cmd),header=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
  if(is.null(d)||nrow(d)<10)return(NULL); names(d)<-c("rsid","b","se"); d[!duplicated(d$rsid)&is.finite(d$b)&is.finite(d$se)&d$se>0,]}
get_narco<-function(rsids){tf<-tempfile(); writeLines(rsids,tf)
  cmd<-sprintf("zcat '%s' | awk 'NR==FNR{a[$1];next}($1 in a){print $1,$6,$7}' '%s' -",NARCO,tf)
  d<-tryCatch(read.table(pipe(cmd),header=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
  unlink(tf); if(is.null(d)||nrow(d)<10)return(NULL); names(d)<-c("rsid","b","se"); d[!duplicated(d$rsid)&is.finite(d$b)&is.finite(d$se)&d$se>0,]}

hits<-read.csv(LAVA,stringsAsFactors=FALSE); hits<-hits[hits$p<1e-3,]
cat(sprintf("narco LAVA hits p<1e-3: %d\n",nrow(hits))); res<-list()
for(k in 1:nrow(hits)){h<-hits[k,]; anc<-ifelse(grepl("^rls",h$pair),"rls_akcimen","pd_nalls")
  a<-get_anchor(anc,h$chr,h$start,h$stop); if(is.null(a)){cat("  no anchor SNPs:",h$pair,h$loc,"\n");next}
  b<-get_narco(a$rsid); if(is.null(b)){cat("  no narco SNPs:",h$pair,h$loc,"\n");next}
  cm<-intersect(a$rsid,b$rsid); if(length(cm)<20){cat("  <20 shared:",h$pair,h$loc,length(cm),"\n");next}
  A<-a[match(cm,a$rsid),]; B<-b[match(cm,b$rsid),]
  d1<-list(beta=A$b,varbeta=A$se^2,snp=cm,type="cc",s=acfg[[anc]]$s,N=acfg[[anc]]$N)
  d2<-list(beta=B$b,varbeta=B$se^2,snp=cm,type="cc",s=NARCO_S,N=NARCO_N)
  r<-tryCatch(suppressWarnings(coloc.abf(d1,d2)),error=function(e){cat("  coloc err:",conditionMessage(e),"\n");NULL}); if(is.null(r))next
  s<-r$summary
  res[[length(res)+1]]<-data.frame(pair=h$pair,region=sprintf("chr%s:%.1f-%.1fMb",h$chr,h$start/1e6,h$stop/1e6),
    lava_rho=round(h$rho,3),lava_p=signif(h$p,2),nSNP=length(cm),
    PP3=round(s["PP.H3.abf"],3),PP4=round(s["PP.H4.abf"],3),
    verdict=ifelse(s["PP.H4.abf"]>0.9,"COLOC",ifelse(s["PP.H4.abf"]>0.8,"suggestive",ifelse(s["PP.H4.abf"]>0.5,"weak","diffuse"))))}
if(length(res)){R<-do.call(rbind,res); R<-R[order(-R$PP4),]; write.table(R,OUT,sep="\t",row.names=FALSE,quote=FALSE); print(R,row.names=FALSE)} else cat("no colocs computed\n")
cat("\n=== NARCO_COLOC_DONE ===\n")
