#!/usr/bin/env Rscript
# Colocalization screen (coloc.abf) of every local genetic correlation reaching the LAVA screening
# threshold: P < 1e-3 in results/lava_test_family.tsv (74 tests; written by lava_test_family.R).
#   - Each GWAS is read over the LAVA block in its own genome build: the GRCh37 block for hg19 files and
#     its GRCh38 coordinates (region_hg38 in results/coloc_screen_full_supplementary.csv) for hg38 files.
#   - Full summary statistics are read for every trait (no HapMap3 restriction).
#   - Narcolepsy (Ollila 2023) is matched to the other trait's block variants by rsID.
#   - Obstructive sleep apnea (GCST90693193, hg38) carries no rsIDs and is matched by position and
#     allele pair.
#   - Three colocalization priors: p12 = 1e-6, 1e-5 (coloc default) and 5e-5.
# For each pair the number of variants read for each trait is recorded (nA, nB); density_flag marks pairs
# in which one trait has fewer than 30% of the other trait's variants in the block.
# Output: results/coloc_screen_blocks.tsv under the analysis root (input to build_table_s3.py).
suppressMessages(library(coloc)); options(scipen = 999)
ROOT <- "/path/to/analysis"; G <- file.path(ROOT, "data/gwas")
REPO <- Sys.getenv("REPO", ".")
OUT <- file.path(ROOT, "results/coloc_screen_blocks.tsv")
CACHE <- file.path(ROOT, "results/coloc_cache"); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
TESTS <- file.path(REPO, "results/lava_test_family.tsv")
BLOCKS <- file.path(REPO, "results/coloc_screen_full_supplementary.csv")   # block coordinates (hg19 and hg38)

## ---- 0. schizophrenia: strip the VCF-style header of the PGC3 release once ----
SCZ <- file.path(G, "psych/scz_eur_clean.tsv.gz")
if (!file.exists(SCZ)) { cat("writing the cleaned schizophrenia file...\n")
  system(sprintf("zcat %s | grep -vE '^##' | sed '1s/^#//' | gzip > %s",
                 shQuote(file.path(G, "psych/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz")), shQuote(SCZ))) }

## ---- 1. tests: P < 1e-3 from the unique test family; hg38 block coordinates ----
U <- read.delim(TESTS, stringsAsFactors = FALSE); H <- U[U$p < 1e-3, ]
S3 <- read.csv(BLOCKS, stringsAsFactors = FALSE)
lab <- c(rls_akcimen="RLS", pd_nalls="PD", anx="anxiety", insomnia="insomnia", mdd="MDD", osa="OSA", ptsd="PTSD", rbd="RBD",
         sleepdur_eur_short="short sleep", sleepdur_eur_long="long sleep", bip="bipolar", chronotype="chronotype",
         daytime_sleepiness="daytime sleepiness", scz="schizophrenia", et="essential tremor", tourette="Tourette", narcolepsy_ollila="narcolepsy")
H$k <- paste(lab[H$phen1], lab[H$phen2], sprintf("chr%s:%s-%s", H$chr, H$start, H$stop))
S3$k <- paste(S3$Trait_1, S3$Trait_2, S3$Region_hg19)
stopifnot(nrow(H) == 74, all(H$k %in% S3$k))
h38 <- do.call(rbind, strsplit(S3$region_hg38[match(H$k, S3$k)], "[:-]"))
H$lo38 <- as.numeric(h38[, 2]); H$hi38 <- as.numeric(h38[, 3]); H$tag <- sprintf("h%02d", seq_len(nrow(H)))

## ---- 2. file definitions ----
# b = genome build; cc/pc = chromosome/position columns; ic/ec/oc = rsID / effect allele / other allele
# (ic = 0: no rsID column); c4/c5 = statistic columns, read according to fl; fc = effect-allele frequency
# (0 = none); ty = cc (case-control, s = case fraction) or quant (N)
F <- list(
 rls_akcimen = list(f="rls/akcimen2026_EUR_RLS.txt", rd="cat", b=38, cc=2, pc=3, ic=10, ec=4, oc=5, c4=6, c5=7, fc=0, fl="beta", ty="cc", s=0.0547),
 pd_nalls    = list(f="pd/nalls2019_PD_eur.h.tsv.gz", rd="zcat", b=38, cc=1, pc=2, ic=9, ec=3, oc=4, c4=5, c5=6, fc=7, fl="beta", ty="cc", s=0.0697),
 insomnia    = list(f="sleep/insomnia_watanabe/Watanabe_2022/insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz", rd="zcat", b=19, cc=2, pc=3, ic=15, ec=4, oc=13, c4=7, c5=8, fc=20, fl="or", ty="cc", s=0.283),
 mdd         = list(f="daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz", rd="zcat", b=19, cc=1, pc=3, ic=2, ec=4, oc=5, c4=9, c5=10, fc=7, fl="or", ty="cc", s=0.218),
 anx         = list(f="psych/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz", rd="zcat", b=19, cc=1, pc=3, ic=2, ec=4, oc=5, c4=9, c5=10, fc=7, fl="or", ty="cc", s=0.143),
 ptsd        = list(f="psych/eur_ptsd_pcs_v4_aug3_2021.vcf.gz", rd="zcat", b=19, cc=1, pc=3, ic=2, ec=4, oc=5, c4=8, c5=6, fc=7, fl="z", ty="cc", s=0.5),
 rbd         = list(f="sleep/rbd_krohn/GCST90204200_buildGRCh37.tsv", rd="cat", b=19, cc=1, pc=2, ic=3, ec=4, oc=5, c4=6, c5=7, fc=8, fl="beta", ty="cc", s=0.112),
 scz         = list(f="psych/scz_eur_clean.tsv.gz", rd="zcat", b=19, cc=1, pc=3, ic=2, ec=4, oc=5, c4=9, c5=10, fc=7, fl="beta", ty="cc", s=0.409),
 bip         = list(f="psych/bip2024_eur_noUKB_no23andMe.gz", rd="zcat", b=19, cc=2, pc=3, ic=1, ec=4, oc=5, c4=7, c5=8, fc=18, fl="or", ty="cc", s=0.072),
 chronotype  = list(f="sleep/chronotype_jones/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz", rd="zcat", b=19, cc=2, pc=3, ic=1, ec=4, oc=5, c4=8, c5=9, fc=6, fl="beta", ty="quant", N=449734),
 daytime_sleepiness = list(f="sleep/daytime_saxena/Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt", rd="cat", b=19, cc=2, pc=3, ic=1, ec=4, oc=5, c4=8, c5=9, fc=6, fl="beta", ty="quant", N=452071),
 sleepdur_eur_short = list(f="sleep/sleepdur_gelernter/eur_meta_short_dbgap.txt.gz", rd="zcat", b=19, cc=2, pc=3, ic=1, ec=4, oc=5, c4=8, c5=9, fc=6, fl="beta_p", ty="cc", s=0.27),
 sleepdur_eur_long  = list(f="sleep/sleepdur_gelernter/eur_meta_long_dbgap.txt.gz",  rd="zcat", b=19, cc=2, pc=3, ic=1, ec=4, oc=5, c4=8, c5=9, fc=6, fl="beta_p", ty="cc", s=0.27),
 osa         = list(f="sleep/osa/GCST90693193.tsv", rd="cat", b=38, cc=1, pc=2, ic=0, ec=3, oc=4, c4=5, c5=6, fc=7, fl="beta", ty="cc", s=0.15))
NARCO <- file.path(G, "sleep/narcolepsy_ollila/CHR_Whites_info0.7.noHLA.noSig.hetP1e-5.LDSC.N.txt.gz")

## ---- 3. one pass per file: extract every block needed from that file ----
extract <- function(key, regs) {   # regs: data.frame(tag, chr, lo, hi); output columns: tag id ea oa c4 c5 eaf pos
  out <- file.path(CACHE, paste0(key, ".tsv")); if (file.exists(out)) return(out)
  cfg <- F[[key]]; stopifnot(!is.null(cfg))
  rf <- tempfile(); write.table(regs, rf, quote = FALSE, row.names = FALSE, col.names = FALSE)
  awk <- sprintf("awk -v cc=%d -v pc=%d -v ic=%d -v ec=%d -v oc=%d -v c4=%d -v c5=%d -v fc=%d 'NR==FNR{n[$2]++; T[$2,n[$2]]=$1; L[$2,n[$2]]=$3; H[$2,n[$2]]=$4; next} {c=$cc; if(!(c in n)) next; p=$pc+0; for(i=1;i<=n[c];i++) if(p>=L[c,i]&&p<=H[c,i]) print T[c,i], (ic>0?$ic:c\":\"p), $ec, $oc, $c4, $c5, (fc>0?$fc:\"NA\"), p}' %s -",
                 cfg$cc, cfg$pc, cfg$ic, cfg$ec, cfg$oc, cfg$c4, cfg$c5, cfg$fc, rf)
  t0 <- Sys.time(); system(sprintf("%s %s 2>/dev/null | %s > %s", cfg$rd, shQuote(file.path(G, cfg$f)), awk, shQuote(out)))
  cat(sprintf("  extracted %-22s %6d rows  (%.0f s)\n", key, as.integer(system(sprintf("wc -l < %s", shQuote(out)), intern = TRUE)), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out }

prep <- function(d, cfg) {   # statistics -> beta / varbeta / MAF
  suppressWarnings({ d$c4 <- as.numeric(d$c4); d$c5 <- as.numeric(d$c5); d$eaf <- as.numeric(d$eaf) })
  qn <- function(p) qnorm(1 - pmin(pmax(p, 1e-300), 0.9999999) / 2)
  if (cfg$fl == "z") { fq <- pmin(d$c5, 1 - d$c5); neff <- d$eaf; vb <- 1 / (2 * fq * (1 - fq) * (neff + d$c4^2)); b <- d$c4 * sqrt(vb); maf <- fq
  } else { if (cfg$fl == "beta") { b <- d$c4; se <- d$c5 } else if (cfg$fl == "or") { b <- log(d$c4); se <- d$c5 } else { b <- d$c4; se <- abs(d$c4) / qn(d$c5) }
           vb <- se^2; maf <- pmin(d$eaf, 1 - d$eaf) }
  d$b <- as.numeric(b); d$vb <- as.numeric(vb); d$maf <- as.numeric(maf)
  d[is.finite(d$b) & is.finite(d$vb) & d$vb > 0, ] }
mkD <- function(d, cfg) { x <- list(beta = d$b, varbeta = d$vb, snp = d$id, MAF = d$maf)
  if (cfg$ty == "quant") { x$type <- "quant"; x$N <- cfg$N } else { x$type <- "cc"; x$s <- cfg$s }; x }
readc <- function(path) { nm <- c("tag", "id", "ea", "oa", "c4", "c5", "eaf", "pos")
  if (file.size(path) == 0) { d <- as.data.frame(setNames(replicate(8, character(0), simplify = FALSE), nm)); d$pos <- numeric(0); return(d) }
  d <- read.table(path, header = FALSE, stringsAsFactors = FALSE, colClasses = "character", fill = TRUE, quote = "", comment.char = "")
  names(d) <- nm; d$pos <- as.numeric(d$pos); d }
# order: rsID filter, then first occurrence of each rsID kept, then statistic conversion
clean <- function(d, cfg, rsonly = TRUE) { if (rsonly) d <- d[grepl("^rs", d$id), ]; d <- d[!duplicated(d$id), ]; prep(d, cfg) }

## ---- 4. block list per file ----
regs <- list()
for (i in seq_len(nrow(H))) { h <- H[i, ]
  for (tr in c(h$phen1, h$phen2)) { if (tr == "narcolepsy_ollila") next
    cfg <- F[[tr]]
    if (cfg$b == 38) { lo <- h$lo38; hi <- h$hi38 } else { lo <- h$start; hi <- h$stop }
    regs[[length(regs) + 1]] <- data.frame(key = tr, tag = h$tag, chr = h$chr, lo = lo, hi = hi) } }
RG <- do.call(rbind, regs)

cat("\n===== extraction =====\n")
paths <- list()
for (key in unique(RG$key)) paths[[key]] <- extract(key, RG[RG$key == key, c("tag", "chr", "lo", "hi")])

## ---- 5. coloc.abf per test ----
res <- list()
for (i in seq_len(nrow(H))) { h <- H[i, ]; c1 <- F[[h$phen1]]
  a <- clean(subset(readc(paths[[h$phen1]]), tag == h$tag), c1)
  if (h$phen2 == "narcolepsy_ollila") {   # narcolepsy: look up the anchor trait's block rsIDs
    tf <- tempfile(); writeLines(a$id, tf)
    nb <- read.table(pipe(sprintf("zcat %s | awk 'NR==FNR{x[$1];next}($1 in x){print $1,$6,$7}' %s -", shQuote(NARCO), tf)), col.names = c("id", "b", "se"), stringsAsFactors = FALSE)
    nb <- nb[!duplicated(nb$id) & is.finite(nb$b) & is.finite(nb$se) & nb$se > 0, ]
    cm <- intersect(a$id, nb$id); A <- a[match(cm, a$id), ]; B <- nb[match(cm, nb$id), ]
    d1 <- list(beta = A$b, varbeta = A$vb, snp = cm, type = "cc", s = c1$s, N = ifelse(h$phen1 == "rls_akcimen", 676175, 482730))
    d2 <- list(beta = B$b, varbeta = B$se^2, snp = cm, type = "cc", s = 0.05, N = 61727); nA <- nrow(a); nB <- nrow(nb)
  } else {
    c2 <- F[[h$phen2]]
    b <- subset(readc(paths[[h$phen2]]), tag == h$tag)
    if (h$phen2 == "osa") {   # OSA: no rsIDs; match to the anchor trait by hg38 position + allele pair
      a$pk <- paste(a$pos, pmin(toupper(a$ea), toupper(a$oa)), pmax(toupper(a$ea), toupper(a$oa)))
      b$pk <- paste(b$pos, pmin(toupper(b$ea), toupper(b$oa)), pmax(toupper(b$ea), toupper(b$oa)))
      b <- b[!duplicated(b$pk), ]; nOSA <- nrow(b); b$id <- a$id[match(b$pk, a$pk)]; b <- b[!is.na(b$id), ]
      cat(sprintf("  OSA position match %s: %d / %d\n", h$k, nrow(b), nOSA)) }
    b <- clean(b, c2)
    nA <- nrow(a); nB <- nrow(b)
    cm <- intersect(a$id, b$id); A <- a[match(cm, a$id), ]; B <- b[match(cm, b$id), ]
    maf <- B$maf; maf[!is.finite(maf)] <- A$maf[!is.finite(maf)]; ok <- is.finite(maf) & maf > 0.01
    A <- A[ok, ]; B <- B[ok, ]; A$maf <- maf[ok]; B$maf <- maf[ok]; cm <- A$id
    d1 <- mkD(A, c1); d2 <- mkD(B, c2) }
  if (length(cm) < 20) { cat("  fewer than 20 shared variants:", h$k, "\n"); next }
  row <- data.frame(k = h$k, nA = nA, nB = nB, nSNP = length(cm))
  for (p12 in c(1e-6, 1e-5, 5e-5)) { r <- suppressWarnings(coloc.abf(d1, d2, p12 = p12)); s <- r$summary
    tag <- if (p12 == 1e-5) "" else if (p12 == 1e-6) "_strict" else "_lenient"
    row[[paste0("PP3", tag)]] <- round(as.numeric(s["PP.H3.abf"]), 3); row[[paste0("PP4", tag)]] <- round(as.numeric(s["PP.H4.abf"]), 3)
    if (p12 == 1e-5) row$lead <- r$results$snp[which.max(r$results$SNP.PP.H4)] }
  row$density_flag <- ifelse(min(nA, nB) / max(nA, nB) < 0.3, "flag", "")
  res[[length(res) + 1]] <- row }
R <- do.call(rbind, res)
write.table(R, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cut3 <- cut(R$PP4, c(-1, 0.8, 0.9, 2), c("none", "suggestive", "COLOC"))
cat(sprintf("\n%d tests screened | PP.H4 > 0.90: %d | 0.80-0.90: %d -> %s\n", nrow(R), sum(cut3 == "COLOC"), sum(cut3 == "suggestive"), OUT))
cat("=== COLOC_SCREEN_FULL_DONE ===\n")
