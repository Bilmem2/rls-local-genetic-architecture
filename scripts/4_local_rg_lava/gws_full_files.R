#!/usr/bin/env Rscript
# Genome-wide-significance class of each screened block, re-evaluated from the full summary statistics for
# short sleep, long sleep and obstructive sleep apnea (known_vs_novel_blocks.R reads these three traits
# from their HapMap3-restricted munge inputs). Other rows are unchanged.
#   short / long sleep (hg19, with rsIDs): block membership by rsID against 1000 Genomes EUR (hg19), as in
#                                          known_vs_novel_blocks.R
#   OSA (hg38, no rsIDs): block membership by position within the block's hg38 coordinates
# Input and output: results/known_vs_novel_blocks.csv (GWS column of Supplementary Table S3).
suppressMessages(library(data.table))
ROOT <- "/path/to/analysis"; G <- file.path(ROOT, "data/gwas")
REPO <- Sys.getenv("REPO", "."); KF <- file.path(REPO, "results/known_vs_novel_blocks.csv")
K <- fread(KF)
S3 <- fread(file.path(REPO, "results/coloc_screen_full_supplementary.csv"))
K <- merge(K, S3[, .(Trait_1, Trait_2, region_hg38, Region_hg19)], by = c("Trait_1", "Trait_2", "region_hg38"), sort = FALSE)
stopifnot(nrow(K) == 74)
bim <- fread(file.path(ROOT, "data/reference/g1000_eur.bim"), header = FALSE, select = c(1, 2, 4), col.names = c("chr", "rsid", "pos"))
FULL <- list("short sleep" = "sleep/sleepdur_gelernter/eur_meta_short_dbgap.txt.gz", "long sleep" = "sleep/sleepdur_gelernter/eur_meta_long_dbgap.txt.gz")
minp_rs <- function(f, chr, lo, hi) {   # sleep-duration files: rsid CHR BP A1 A2 Freq1 N Effect P
  d <- fread(cmd = sprintf("zcat %s | awk 'NR>1&&$2==%d&&$3>=%.0f&&$3<=%.0f&&$1~/^rs/{print $1,$9}'", shQuote(file.path(G, f)), chr, lo, hi), header = FALSE, col.names = c("rsid", "p"))
  d <- d[rsid %in% bim[chr == chr & pos >= lo & pos <= hi, rsid]]
  if (!nrow(d)) return(list(NA, NA)); i <- which.min(d$p); list(d$p[i], d$rsid[i]) }
minp_osa <- function(chr, lo, hi) {     # OSA file: chromosome base_pair_location EA OA beta se eaf p (hg38)
  d <- fread(cmd = sprintf("awk 'NR>1&&$1==%d&&$2>=%.0f&&$2<=%.0f{print $2,$8}' %s", chr, lo, hi, shQuote(file.path(G, "sleep/osa/GCST90693193.tsv"))), header = FALSE, col.names = c("pos", "p"))
  if (!nrow(d)) return(list(NA, NA)); i <- which.min(d$p); list(d$p[i], paste0(chr, ":", d$pos[i])) }
for (i in which(K$Trait_2 %in% c("short sleep", "long sleep", "OSA"))) {
  r19 <- as.numeric(strsplit(sub("^chr", "", K$Region_hg19[i]), "[:-]")[[1]]); r38 <- as.numeric(strsplit(K$region_hg38[i], "[:-]")[[1]])
  v <- if (K$Trait_2[i] == "OSA") minp_osa(r38[1], r38[2], r38[3]) else minp_rs(FULL[[K$Trait_2[i]]], r19[1], r19[2], r19[3])
  K$minp_t2[i] <- v[[1]]; K$lead_t2[i] <- v[[2]] }
gws1 <- K$minp_t1 < 5e-8; gws2 <- K$minp_t2 < 5e-8
K[, class := fifelse(Trait_2 == "narcolepsy", class, fifelse(gws1 & gws2, "both", fifelse(gws1, "trait 1 only", fifelse(gws2, "trait 2 only", "neither"))))]
K <- K[, .(Trait_1, Trait_2, Region_hg19, region_hg38, minp_t1, lead_t1, minp_t2, lead_t2, class)]
cat("Classes:\n"); print(table(K$class))
fwrite(K, KF); cat("written:", KF, "\n")
