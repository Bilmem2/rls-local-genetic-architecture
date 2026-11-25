#!/usr/bin/env Rscript
# For each of the 74 screened local genetic correlations (Supplementary Table S3), ask whether
# either trait's own GWAS carries a genome-wide-significant variant (P < 5e-8) inside the LAVA block.
#   both     -> known locus for both traits
#   t1 / t2  -> known locus for one trait only
#   neither  -> no genome-wide-significant variant for either trait in the GWAS analysed
# Block membership is by rsID against the 1000 Genomes EUR reference used by LAVA (hg19), so the
# hg38 files (RLS Akcimen, PD Nalls) and the hg19 files are treated identically.
# Narcolepsy (Ollila 2023) is not assessable: its public release removes every variant with P < 5e-8.
# The short-sleep, long-sleep and OSA rows are then re-evaluated from the full files by gws_full_files.R.
suppressMessages(library(data.table))
ROOT <- "/path/to/analysis"; G <- file.path(ROOT, "data/gwas")
REPO <- Sys.getenv("REPO", ".")  # repository root; run_all.sh runs from there
HITS <- Sys.getenv("HITS", file.path(REPO, "results/coloc_screen_full_supplementary.csv"))
BIM  <- file.path(ROOT, "data/reference/g1000_eur.bim")
OUT  <- Sys.getenv("OUT", file.path(REPO, "results/known_vs_novel_blocks.csv"))
GWS  <- 5e-8; PAD <- 3e6

# trait configs: file, reader, type (cc/quant), s or N; awk prints rsid ea oa c4 c5 eaf
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
# if the cleaned schizophrenia file is absent, read the PGC3 original (same column layout)
if (!file.exists(file.path(G, cf$scz$f))) cf$scz$f <- "psych/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"

key <- c("RLS" = "rls_akcimen", "PD" = "pd_nalls", "insomnia" = "insomnia", "MDD" = "mdd",
         "anxiety" = "anx", "PTSD" = "ptsd", "OSA" = "osa", "short sleep" = "sleepdur_eur_short",
         "long sleep" = "sleepdur_eur_long", "RBD" = "rbd", "schizophrenia" = "scz",
         "bipolar" = "bip", "chronotype" = "chronotype", "daytime sleepiness" = "daytime_sleepiness",
         "narcolepsy" = NA)

H <- fread(HITS)
rg <- tstrsplit(sub("^chr", "", H$Region_hg19), "[:-]")
H[, `:=`(chr = as.integer(rg[[1]]), start = as.numeric(rg[[2]]), stop = as.numeric(rg[[3]]))]
H[, blk := paste0(chr, ":", start, "-", stop)]
B <- unique(H[, .(blk, chr, start, stop)])
cat(sprintf("rows %d | unique blocks %d\n", nrow(H), nrow(B)))

# rsIDs of the reference panel inside each block (hg19)
bim <- fread(BIM, header = FALSE, select = c(1, 2, 4), col.names = c("chr", "rsid", "pos"))
bim <- bim[chr %in% unique(B$chr)]
memb <- rbindlist(lapply(seq_len(nrow(B)), function(i)
  bim[chr == B$chr[i] & pos >= B$start[i] & pos <= B$stop[i], .(rsid, blk = B$blk[i])]))
rm(bim); invisible(gc())
cat(sprintf("reference rsIDs inside blocks: %d\n", nrow(memb)))

zfrom <- function(d, fl) {
  qn <- function(p) qnorm(1 - pmin(pmax(p, 1e-300), 0.9999999) / 2)
  suppressWarnings({ c4 <- as.numeric(d$c4); c5 <- as.numeric(d$c5) })
  if (fl == "z") return(c4)
  if (fl == "beta") return(c4 / c5)
  if (fl == "or") return(log(c4) / c5)
  if (fl %in% c("beta_p", "or_p")) return(qn(c5) * sign(if (fl == "or_p") log(c4) else c4))
  stop("unknown flag ", fl)
}

minp <- list()
for (lab in setdiff(names(key), "narcolepsy")) {
  tr <- key[[lab]]; cc <- cf[[tr]]
  rows <- H[Trait_1 == lab | Trait_2 == lab]
  if (!nrow(rows)) next
  W <- unique(rows[, .(chr, lo = pmax(0, start - PAD), hi = stop + PAD)])
  m <- regmatches(cc$awk, regexec("\\$(\\d+)==c&&\\$(\\d+)>=lo&&\\$(\\d+)<=hi", cc$awk))[[1]]
  if (length(m) != 4) stop("cannot parse awk for ", tr)
  win <- paste(sprintf("($%s==%d&&$%s>=%.0f&&$%s<=%.0f)", m[2], W$chr, m[3], W$lo, m[3], W$hi), collapse = "||")
  awk <- sub(m[1], paste0("(", win, ")"), cc$awk, fixed = TRUE)
  cmd <- sprintf("%s %s 2>/dev/null | awk '%s'", cc$rd, shQuote(file.path(G, cc$f)), awk)
  d <- tryCatch(fread(cmd = cmd, header = FALSE, fill = TRUE, col.names = c("rsid", "ea", "oa", "c4", "c5", "eaf")),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) { cat(sprintf("%-20s NO DATA\n", lab)); next }
  d <- unique(d, by = "rsid")
  d[, p := 2 * pnorm(-abs(zfrom(.SD, cc$fl)))]
  d <- merge(d[is.finite(p), .(rsid, p)], memb, by = "rsid")
  mp <- d[, .(minp = min(p), lead = rsid[which.min(p)], n = .N), by = blk]
  mp[, trait := lab]
  minp[[lab]] <- mp
  cat(sprintf("%-20s %7d variants in blocks | %d blocks\n", lab, nrow(d), nrow(mp)))
}
MP <- rbindlist(minp)

get <- function(t, b, col) { v <- MP[trait == t & blk == b][[col]]; if (length(v)) v else NA }
H[, `:=`(minp_t1 = mapply(get, Trait_1, blk, "minp"), lead_t1 = mapply(get, Trait_1, blk, "lead"),
         minp_t2 = mapply(get, Trait_2, blk, "minp"), lead_t2 = mapply(get, Trait_2, blk, "lead"))]
H[, `:=`(gws_t1 = minp_t1 < GWS, gws_t2 = minp_t2 < GWS)]
H[, class := fifelse(Trait_2 == "narcolepsy",
                     fifelse(gws_t1, "RLS/PD only (narcolepsy n/a)", "neither (narcolepsy n/a)"),
              fifelse(gws_t1 & gws_t2, "both",
              fifelse(gws_t1, "trait 1 only", fifelse(gws_t2, "trait 2 only", "neither"))))]

fwrite(H[, .(Trait_1, Trait_2, region_hg38, LAVA_rho, LAVA_P, verdict, survives_bonferroni,
             minp_t1 = signif(minp_t1, 3), lead_t1, minp_t2 = signif(minp_t2, 3), lead_t2, class)], OUT)
cat("\n== class counts ==\n"); print(H[, .N, by = class][order(-N)])
cat("\n== by verdict ==\n"); print(dcast(H, class ~ verdict, fun.aggregate = length, value.var = "blk"))
cat("\n== neither, strongest first ==\n")
print(H[class == "neither"][order(LAVA_P), .(Trait_1, Trait_2, region_hg38, LAVA_rho, LAVA_P,
                                             minp_t1 = signif(minp_t1, 2), minp_t2 = signif(minp_t2, 2))])
cat("\nwritten:", OUT, "\n")
