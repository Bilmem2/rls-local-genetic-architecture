#!/usr/bin/env Rscript
# =============================================================================
# panel_lava.R — genome-wide LAVA local rg for the RLS/PD x sleep+psych panel.
#   Anchors: rls_akcimen (primary RLS, biggest) + pd_nalls (primary PD).
#   Overlap matrices auto-built from LDSC --rg intercepts (results/ldsc_rg_panel).
#   Panel traits treated continuous (local rho is scale-invariant) -> N from sumstats.
#   Both-trait univariate local-h2 gate (p<0.05/2495) then bivariate local rg.
#   INCREMENTAL write: each pair's hits appended to CSV as it finishes.
# =============================================================================
suppressMessages(library(LAVA))
ROOT <- "/path/to/analysis"
REF  <- file.path(ROOT,"data/reference/g1000_eur")
LOCF <- file.path(ROOT,"data/reference/lava_blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
M    <- file.path(ROOT,"results/munged")
RGDIR<- file.path(ROOT,"results/ldsc_rg_panel")
Ldir <- file.path(ROOT,"results/lava"); dir.create(Ldir, showWarnings=FALSE, recursive=TRUE)
OUT  <- file.path(Ldir,"lava_panel_bivariate.csv")
UNIVOUT <- file.path(Ldir,"lava_panel_univ_gatepass.csv")
loci <- read.loci(LOCF); n <- nrow(loci); N_GENOME <- 2495; ut <- 0.05/N_GENOME
cat(sprintf("PANEL LAVA | loci=%d | both-trait univ gate p<%.2g\n", n, ut)); flush.console()

# binary anchors get case/control counts; every other trait -> continuous (NA NA)
cc <- list(rls_akcimen=c(36993,639182), pd_nalls=c(33674,449056))
info_row <- function(ph){ f<-file.path(M,paste0(ph,".sumstats.gz"))
  if(!is.null(cc[[ph]])) sprintf("%s %d %d %s",ph,cc[[ph]][1],cc[[ph]][2],f) else sprintf("%s NA NA %s",ph,f) }
write_info <- function(file,p1,p2) writeLines(c("phenotype cases controls filename",info_row(p1),info_row(p2)),file)

# parse an LDSC --rg .log: summary table (per-trait h2_int + gcov_int) and anchor's own intercept
read_rg_tab <- function(logf){ if(!file.exists(logf)) return(NULL); L<-readLines(logf,warn=FALSE)
  i<-grep("Summary of Genetic Correlation Results",L); if(!length(i)) return(NULL)
  j<-i[1]+1; k<-j; while(k+1<=length(L) && nchar(trimws(L[k+1]))>0) k<-k+1
  tab<-read.table(text=L[j:k],header=TRUE,stringsAsFactors=FALSE,fill=TRUE)
  tab$trait<-sub("\\.sumstats\\.gz$","",basename(tab$p2)); tab }
anchor_intercept <- function(logf){ L<-readLines(logf,warn=FALSE)
  i<-grep("Heritability of phenotype 1",L)[1]; seg<-L[i:min(i+10,length(L))]
  ln<-grep("Intercept:",seg,value=TRUE)[1]; as.numeric(sub(".*Intercept:\\s*([0-9.]+).*","\\1",ln)) }
build_ovl <- function(file,anchor,trait,tab,aint){ r<-tab[tab$trait==trait,]
  g<-as.numeric(r$gcov_int[1]); ht<-as.numeric(r$h2_int[1])
  m<-matrix(c(aint,g,g,ht),2,dimnames=list(c(anchor,trait),c(anchor,trait)))
  write.table(round(cov2cor(m),5),file,quote=FALSE) }

rls_log<-file.path(RGDIR,"rls_akcimen_vs_panel.log"); pd_log<-file.path(RGDIR,"pd_vs_panel.log")
rls_tab<-read_rg_tab(rls_log); rls_aint<-anchor_intercept(rls_log)
pd_tab <-read_rg_tab(pd_log);  pd_aint <-anchor_intercept(pd_log)
for(ex in c("insomnia","daytime")){ ef<-file.path(RGDIR,paste0("pd_nalls_vs_",ex,".log"))
  et<-read_rg_tab(ef); if(!is.null(et)) pd_tab<-rbind(pd_tab,et) }
cat(sprintf("rls_akcimen intercept=%.4f (%d traits) ; pd_nalls intercept=%.4f (%d traits)\n",
    rls_aint,nrow(rls_tab),pd_aint,nrow(pd_tab))); flush.console()

# priority-ordered pairs: tag, anchor, trait, src(rls|pd)
PAIRS <- list(
  c("rlsak_pd","rls_akcimen","pd_nalls","rls"),
  c("rlsak_insomnia","rls_akcimen","insomnia","rls"),
  c("pd_insomnia","pd_nalls","insomnia","pd"),
  c("pd_anx","pd_nalls","anx","pd"), c("pd_mdd","pd_nalls","mdd","pd"),
  c("pd_ptsd","pd_nalls","ptsd","pd"), c("pd_osa","pd_nalls","osa","pd"),
  c("pd_sleepdur","pd_nalls","sleepdur_eur_short","pd"), c("pd_rbd","pd_nalls","rbd","pd"),
  c("rlsak_anx","rls_akcimen","anx","rls"), c("rlsak_mdd","rls_akcimen","mdd","rls"),
  c("rlsak_ptsd","rls_akcimen","ptsd","rls"), c("rlsak_osa","rls_akcimen","osa","rls"),
  c("rlsak_sleepdur","rls_akcimen","sleepdur_eur_short","rls"), c("rlsak_rbd","rls_akcimen","rbd","rls"),
  c("pd_scz","pd_nalls","scz","pd"), c("pd_bip","pd_nalls","bip","pd"),
  c("pd_chronotype","pd_nalls","chronotype","pd"), c("pd_daytime","pd_nalls","daytime_sleepiness","pd"),
  c("pd_sleepdurlong","pd_nalls","sleepdur_eur_long","pd"), c("pd_narco","pd_nalls","narcolepsy","pd"),
  c("rlsak_scz","rls_akcimen","scz","rls"), c("rlsak_bip","rls_akcimen","bip","rls"),
  c("rlsak_chronotype","rls_akcimen","chronotype","rls"), c("rlsak_daytime","rls_akcimen","daytime_sleepiness","rls"),
  c("rlsak_sleepdurlong","rls_akcimen","sleepdur_eur_long","rls"), c("rlsak_narco","rls_akcimen","narcolepsy","rls"))

first <- TRUE; ufirst <- TRUE
for(P in PAIRS){
  tag<-P[1]; anchor<-P[2]; trait<-P[3]; src<-P[4]
  tab <- if(src=="rls") rls_tab else pd_tab; aint <- if(src=="rls") rls_aint else pd_aint
  if(is.null(tab) || !(trait %in% tab$trait)){ cat("  skip (no rg intercepts):",tag,"\n"); next }
  info<-file.path(Ldir,paste0("info.",tag,".txt")); ovl<-file.path(Ldir,paste0("ovl.",tag,".txt"))
  write_info(info,anchor,trait)
  if(!tryCatch({build_ovl(ovl,anchor,trait,tab,aint);TRUE},error=function(e){cat(" ovl err",tag,conditionMessage(e),"\n");FALSE})) next
  cat(sprintf("\n>>> [%s] %s x %s\n",tag,anchor,trait)); flush.console()
  inp <- tryCatch(process.input(info,ovl,REF,c(anchor,trait)),error=function(e){cat(" input err:",conditionMessage(e),"\n");NULL})
  if(is.null(inp)) next
  res<-list(); np<-0; ng<-0
  for(i in 1:n){
    if(i%%400==0) message(sprintf("   [%s] %d/%d gatepass %d hits %d",tag,i,n,ng,np))
    loc<-tryCatch(process.locus(loci[i,],inp),error=function(e) NULL); if(is.null(loc)) next
    u<-tryCatch(run.univ(loc),error=function(e) NULL); if(is.null(u)) next
    if(any(u$phen==anchor & u$p<ut) && any(u$phen==trait & u$p<ut)){
      ng<-ng+1
      b<-tryCatch(run.bivar(loc,phenos=c(anchor,trait)),error=function(e) NULL)
      if(!is.null(b) && nrow(b)>0){
        b$pair<-tag; b$loc<-as.character(loci$LOC[i]); b$chr<-loc$chr; b$start<-loc$start; b$stop<-loc$stop
        res[[length(res)+1]]<-b; np<-np+1
        cat(sprintf("   [%s] %s chr%s rho=%.3f p=%.2g\n",tag,loci$LOC[i],loc$chr,b$rho,b$p)); flush.console()
      }
    }
  }
  cat(sprintf("   [%s] gate-passed=%d bivariate-hits=%d\n",tag,ng,np))
  if(length(res)>0){ df<-do.call(rbind,res)
    write.table(df,OUT,sep=",",row.names=FALSE,col.names=first,append=!first); first<-FALSE }
}
cat("\n=== PANEL_LAVA_DONE ===\n")
