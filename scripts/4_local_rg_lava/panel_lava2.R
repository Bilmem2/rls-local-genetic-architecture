#!/usr/bin/env Rscript
# panel_lava2.R — parallelizable: PAIR_GROUP env selects a disjoint subset + own CSV.
#   pdsig  = PD x {insomnia,anx,mdd,ptsd,osa,sleepdur,rbd}
#   rlssig = rls_akcimen x {anx,mdd,ptsd,osa,sleepdur,rbd}
#   tail   = both anchors x {scz,bip,chronotype,daytime,sleepdur_long,narcolepsy}
#   (rlsak_pd + rlsak_insomnia already done in lava_panel_bivariate.csv)
suppressMessages(library(LAVA))
GRP <- Sys.getenv("PAIR_GROUP","pdsig")
ROOT <- "/path/to/analysis"
REF  <- file.path(ROOT,"data/reference/g1000_eur")
LOCF <- file.path(ROOT,"data/reference/lava_blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
M    <- file.path(ROOT,"results/munged"); RGDIR<- file.path(ROOT,"results/ldsc_rg_panel")
Ldir <- file.path(ROOT,"results/lava")
OUT  <- file.path(Ldir, paste0("lava_panel_",GRP,".csv"))
loci <- read.loci(LOCF); n <- nrow(loci); ut <- 0.05/2495
cat(sprintf("PANEL LAVA [%s] | loci=%d | gate p<%.2g -> %s\n", GRP, n, ut, OUT)); flush.console()

cc <- list(rls_akcimen=c(36993,639182), pd_nalls=c(33674,449056))
info_row <- function(ph){ f<-file.path(M,paste0(ph,".sumstats.gz"))
  if(!is.null(cc[[ph]])) sprintf("%s %d %d %s",ph,cc[[ph]][1],cc[[ph]][2],f) else sprintf("%s NA NA %s",ph,f) }
write_info <- function(file,p1,p2) writeLines(c("phenotype cases controls filename",info_row(p1),info_row(p2)),file)
read_rg_tab <- function(logf){ if(!file.exists(logf)) return(NULL); L<-readLines(logf,warn=FALSE)
  i<-grep("Summary of Genetic Correlation Results",L); if(!length(i)) return(NULL)
  j<-i[1]+1; k<-j; while(k+1<=length(L) && nchar(trimws(L[k+1]))>0) k<-k+1
  tab<-read.table(text=L[j:k],header=TRUE,stringsAsFactors=FALSE,fill=TRUE)
  tab$trait<-sub("\\.sumstats\\.gz$","",basename(tab$p2)); tab }
anchor_intercept <- function(logf){ L<-readLines(logf,warn=FALSE)
  i<-grep("Heritability of phenotype 1",L)[1]; seg<-L[i:min(i+10,length(L))]
  ln<-grep("Intercept:",seg,value=TRUE)[1]; as.numeric(sub(".*Intercept:\\s*([0-9.]+).*","\\1",ln)) }
build_ovl <- function(file,anchor,trait,tab,aint){ r<-tab[tab$trait==trait,]
  m<-matrix(c(aint,as.numeric(r$gcov_int[1]),as.numeric(r$gcov_int[1]),as.numeric(r$h2_int[1])),2,
            dimnames=list(c(anchor,trait),c(anchor,trait))); write.table(round(cov2cor(m),5),file,quote=FALSE) }

rls_tab<-read_rg_tab(file.path(RGDIR,"rls_akcimen_vs_panel.log")); rls_aint<-anchor_intercept(file.path(RGDIR,"rls_akcimen_vs_panel.log"))
pd_tab <-read_rg_tab(file.path(RGDIR,"pd_vs_panel.log")); pd_aint<-anchor_intercept(file.path(RGDIR,"pd_vs_panel.log"))
for(ex in c("insomnia","daytime")){ et<-read_rg_tab(file.path(RGDIR,paste0("pd_nalls_vs_",ex,".log"))); if(!is.null(et)) pd_tab<-rbind(pd_tab,et) }

G <- list(
  pdsig  = list(c("pd_insomnia","pd_nalls","insomnia","pd"),c("pd_anx","pd_nalls","anx","pd"),c("pd_mdd","pd_nalls","mdd","pd"),
                c("pd_ptsd","pd_nalls","ptsd","pd"),c("pd_osa","pd_nalls","osa","pd"),c("pd_sleepdur","pd_nalls","sleepdur_eur_short","pd"),c("pd_rbd","pd_nalls","rbd","pd")),
  rlssig = list(c("rlsak_anx","rls_akcimen","anx","rls"),c("rlsak_mdd","rls_akcimen","mdd","rls"),c("rlsak_ptsd","rls_akcimen","ptsd","rls"),
                c("rlsak_osa","rls_akcimen","osa","rls"),c("rlsak_sleepdur","rls_akcimen","sleepdur_eur_short","rls"),c("rlsak_rbd","rls_akcimen","rbd","rls")),
  tail   = list(c("pd_scz","pd_nalls","scz","pd"),c("pd_bip","pd_nalls","bip","pd"),c("pd_chronotype","pd_nalls","chronotype","pd"),
                c("pd_daytime","pd_nalls","daytime_sleepiness","pd"),c("pd_sleepdurlong","pd_nalls","sleepdur_eur_long","pd"),c("pd_narco","pd_nalls","narcolepsy","pd"),
                c("rlsak_scz","rls_akcimen","scz","rls"),c("rlsak_bip","rls_akcimen","bip","rls"),c("rlsak_chronotype","rls_akcimen","chronotype","rls"),
                c("rlsak_daytime","rls_akcimen","daytime_sleepiness","rls"),c("rlsak_sleepdurlong","rls_akcimen","sleepdur_eur_long","rls"),c("rlsak_narco","rls_akcimen","narcolepsy","rls")))
PAIRS <- G[[GRP]]; first <- TRUE
for(P in PAIRS){ tag<-P[1]; anchor<-P[2]; trait<-P[3]; src<-P[4]
  tab <- if(src=="rls") rls_tab else pd_tab; aint <- if(src=="rls") rls_aint else pd_aint
  if(is.null(tab) || !(trait %in% tab$trait)){ cat("  skip (no rg):",tag,"\n"); next }
  info<-file.path(Ldir,paste0("info.",tag,".txt")); ovl<-file.path(Ldir,paste0("ovl.",tag,".txt"))
  write_info(info,anchor,trait)
  if(!tryCatch({build_ovl(ovl,anchor,trait,tab,aint);TRUE},error=function(e) FALSE)) next
  cat(sprintf("\n>>> [%s] %s x %s\n",tag,anchor,trait)); flush.console()
  inp <- tryCatch(process.input(info,ovl,REF,c(anchor,trait)),error=function(e){cat(" input err\n");NULL}); if(is.null(inp)) next
  res<-list(); np<-0; ng<-0
  for(i in 1:n){ if(i%%400==0) message(sprintf("   [%s] %d/%d gate %d hit %d",tag,i,n,ng,np))
    loc<-tryCatch(process.locus(loci[i,],inp),error=function(e) NULL); if(is.null(loc)) next
    u<-tryCatch(run.univ(loc),error=function(e) NULL); if(is.null(u)) next
    if(any(u$phen==anchor & u$p<ut) && any(u$phen==trait & u$p<ut)){ ng<-ng+1
      b<-tryCatch(run.bivar(loc,phenos=c(anchor,trait)),error=function(e) NULL)
      if(!is.null(b) && nrow(b)>0){ b$pair<-tag; b$loc<-as.character(loci$LOC[i]); b$chr<-loc$chr; b$start<-loc$start; b$stop<-loc$stop
        res[[length(res)+1]]<-b; np<-np+1; cat(sprintf("   [%s] %s chr%s rho=%.3f p=%.2g\n",tag,loci$LOC[i],loc$chr,b$rho,b$p)); flush.console() } } }
  cat(sprintf("   [%s] gate-passed=%d hits=%d\n",tag,ng,np))
  if(length(res)>0){ df<-do.call(rbind,res); write.table(df,OUT,sep=",",row.names=FALSE,col.names=first,append=!first); first<-FALSE } }
cat(sprintf("\n=== PANEL_LAVA_%s_DONE ===\n",GRP))
