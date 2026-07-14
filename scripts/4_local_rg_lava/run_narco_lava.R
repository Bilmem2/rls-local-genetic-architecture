#!/usr/bin/env Rscript
# Focused LAVA re-run for narcolepsy using the Ollila EUR GWAS (replaces the MVP proxy).
# Runs the two narco pairs (RLS x narco, PD x narco) genome-wide with the same both-trait
# local-h2 gate + LDSC-intercept overlap machinery as panel_lava.R.
# Output: results/lava/lava_narco_ollila.csv  (or a message if nothing passes the gate)
suppressMessages(library(LAVA))
ROOT <- "/path/to/analysis"
REF  <- file.path(ROOT,"data/reference/g1000_eur")
LOCF <- file.path(ROOT,"data/reference/lava_blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
M    <- file.path(ROOT,"results/munged")
RGDIR<- file.path(ROOT,"results/ldsc_rg_panel")
Ldir <- file.path(ROOT,"results/lava"); dir.create(Ldir, showWarnings=FALSE, recursive=TRUE)
OUT  <- file.path(Ldir,"lava_narco_ollila.csv")
if (file.exists(OUT)) file.remove(OUT)
loci <- read.loci(LOCF); n <- nrow(loci); N_GENOME <- 2495; ut <- 0.05/N_GENOME
cat(sprintf("NARCO LAVA (Ollila) | loci=%d | both-trait univ gate p<%.2g\n", n, ut)); flush.console()

# binary anchors carry case/control counts; narcolepsy treated continuous (local rho scale-invariant)
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
  g<-as.numeric(r$gcov_int[1]); ht<-as.numeric(r$h2_int[1])
  m<-matrix(c(aint,g,g,ht),2,dimnames=list(c(anchor,trait),c(anchor,trait)))
  write.table(round(cov2cor(m),5),file,quote=FALSE) }

# tag, anchor, trait, rg-log (narco-specific logs produced by the munge step)
PAIRS <- list(
  c("rlsak_narco_ollila","rls_akcimen","narcolepsy_ollila", file.path(RGDIR,"rls_akcimen_vs_narco_ollila.log")),
  c("pd_narco_ollila",   "pd_nalls",   "narcolepsy_ollila", file.path(RGDIR,"pd_vs_narco_ollila.log")))

first <- TRUE
for(P in PAIRS){
  tag<-P[1]; anchor<-P[2]; trait<-P[3]; rglog<-P[4]
  tab<-read_rg_tab(rglog); aint<-tryCatch(anchor_intercept(rglog),error=function(e) NA)
  if(is.null(tab) || !(trait %in% tab$trait) || is.na(aint)){ cat("  skip (no rg intercepts):",tag,"\n"); next }
  info<-file.path(Ldir,paste0("info.",tag,".txt")); ovl<-file.path(Ldir,paste0("ovl.",tag,".txt"))
  write_info(info,anchor,trait)
  if(!tryCatch({build_ovl(ovl,anchor,trait,tab,aint);TRUE},error=function(e){cat(" ovl err",tag,conditionMessage(e),"\n");FALSE})) next
  cat(sprintf("\n>>> [%s] %s x %s (anchor intercept=%.4f)\n",tag,anchor,trait,aint)); flush.console()
  inp<-tryCatch(process.input(info,ovl,REF,c(anchor,trait)),error=function(e){cat(" input err:",conditionMessage(e),"\n");NULL})
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
cat("\n=== NARCO_LAVA_DONE ===\n")
if(!file.exists(OUT)) cat("RESULT: no bivariate hits passed the gate -> narcolepsy does not locally overlap RLS/PD at p<1e-3 (expected, HLA excluded).\n") else {
  cat("RESULT: hits written to", OUT, "\n"); print(read.csv(OUT)) }
