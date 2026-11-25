#!/usr/bin/env Rscript
# Robustness of the eight colocalizing pairs of the screen (seven with PP.H4 > 0.90, one at 0.80-0.90):
#   (1) coloc.abf over windows of +-100 kb, +-250 kb, +-500 kb and +-1 Mb around the lead variant, each
#       under p12 = 1e-6, 1e-5 and 5e-5 (Supplementary Table S4, panels A and B);
#   (2) coloc-SuSiE over +-250 kb around the lead variant, with LD from the 1000 Genomes Phase 3 European
#       samples distributed with LDSC (n = 489; PLINK v1.9), runsusie max_iter 300 (panel C).
# File definitions, readers and statistic conversions are those of coloc_screen_full.R; windows are taken
# around the lead variant's position in each file's own genome build.
# Outputs (results/): coloc_window_prior_sensitivity.tsv, coloc_susie_summary.tsv, susie_<pair>.csv
suppressMessages({ library(coloc); library(susieR) }); options(scipen = 999)
ROOT <- "/path/to/analysis"; G <- file.path(ROOT, "data/gwas")
REPO <- Sys.getenv("REPO", "."); RES <- file.path(REPO, "results")
CACHE <- file.path(ROOT, "results/coloc_sens_cache"); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
WD <- file.path(ROOT, "results/susie_ld"); dir.create(WD, showWarnings = FALSE, recursive = TRUE)
PLINK <- Sys.getenv("PLINK", "plink"); BDIR <- file.path(ROOT, "data/ldsc_ref/1000G_EUR_Phase3_plink")

## ---- file definitions (as in coloc_screen_full.R; N is the sample size passed to SuSiE) ----
F <- list(
 rls_akcimen = list(f="rls/akcimen2026_EUR_RLS.txt", rd="cat", b=38, cc=2, pc=3, ic=10, ec=4, oc=5, c4=6, c5=7, fc=0, fl="beta", ty="cc", s=0.0547, N=676175),
 pd_nalls    = list(f="pd/nalls2019_PD_eur.h.tsv.gz", rd="zcat", b=38, cc=1, pc=2, ic=9, ec=3, oc=4, c4=5, c5=6, fc=7, fl="beta", ty="cc", s=0.0697, N=482730),
 insomnia    = list(f="sleep/insomnia_watanabe/Watanabe_2022/insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz", rd="zcat", b=19, cc=2, pc=3, ic=15, ec=4, oc=13, c4=7, c5=8, fc=20, fl="or", ty="cc", s=0.283, N=386988),
 mdd         = list(f="daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz", rd="zcat", b=19, cc=1, pc=3, ic=2, ec=4, oc=5, c4=9, c5=10, fc=7, fl="or", ty="cc", s=0.218, N=1639572),
 rbd         = list(f="sleep/rbd_krohn/GCST90204200_buildGRCh37.tsv", rd="cat", b=19, cc=1, pc=2, ic=3, ec=4, oc=5, c4=6, c5=7, fc=8, fl="beta", ty="cc", s=0.112, N=9447),
 chronotype  = list(f="sleep/chronotype_jones/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz", rd="zcat", b=19, cc=2, pc=3, ic=1, ec=4, oc=5, c4=8, c5=9, fc=6, fl="beta", ty="quant", N=449734),
 daytime_sleepiness = list(f="sleep/daytime_saxena/Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt", rd="cat", b=19, cc=2, pc=3, ic=1, ec=4, oc=5, c4=8, c5=9, fc=6, fl="beta", ty="quant", N=452071),
 sleepdur_eur_short = list(f="sleep/sleepdur_gelernter/eur_meta_short_dbgap.txt.gz", rd="zcat", b=19, cc=2, pc=3, ic=1, ec=4, oc=5, c4=8, c5=9, fc=6, fl="beta_p", ty="cc", s=0.27, N=445966),
 osa         = list(f="sleep/osa/GCST90693193.tsv", rd="cat", b=38, cc=1, pc=2, ic=0, ec=3, oc=4, c4=5, c5=6, fc=7, fl="beta", ty="cc", s=0.15, N=980757))

PAIRS <- data.frame(
 id   = c("meis1_insomnia", "meis1_chronotype", "meis1_short_sleep", "skor2_daytime_sleepiness", "sema6d_mdd", "tox3_pd", "ccdc148_pkp4_osa", "gba1_rbd"),
 t1   = c("rls_akcimen","rls_akcimen","rls_akcimen","rls_akcimen","rls_akcimen","rls_akcimen","rls_akcimen","pd_nalls"),
 t2   = c("insomnia","chronotype","sleepdur_eur_short","daytime_sleepiness","mdd","pd_nalls","osa","rbd"),
 chr  = c(2, 2, 2, 18, 15, 16, 2, 1),
 lead = c("rs113851554","rs113851554","rs113851554","rs9965170","rs586226","rs3104769","rs77571593","rs35749011"),
 p19  = c(66750564, 66750564, 66750564, 44788274, 48018198, 52627231, 159300402, 155135036),
 p38  = c(66523432, 66523432, 66523432, 47261903, 47726001, 52593319, 158443890, 155162560), stringsAsFactors = FALSE)
WMAX <- 1e6; WINS <- c(1e5, 2.5e5, 5e5, 1e6); PRI <- c(strict = 1e-6, default = 1e-5, lenient = 5e-5)

## ---- readers (as in coloc_screen_full.R) ----
extract <- function(key, regs) {
  out <- file.path(CACHE, paste0(key, ".tsv")); if (file.exists(out)) return(out); cfg <- F[[key]]
  rf <- tempfile(); write.table(regs, rf, quote = FALSE, row.names = FALSE, col.names = FALSE)
  awk <- sprintf("awk -v cc=%d -v pc=%d -v ic=%d -v ec=%d -v oc=%d -v c4=%d -v c5=%d -v fc=%d 'NR==FNR{n[$2]++; T[$2,n[$2]]=$1; L[$2,n[$2]]=$3; H[$2,n[$2]]=$4; next} {c=$cc; if(!(c in n)) next; p=$pc+0; for(i=1;i<=n[c];i++) if(p>=L[c,i]&&p<=H[c,i]) print T[c,i], (ic>0?$ic:c\":\"p), $ec, $oc, $c4, $c5, (fc>0?$fc:\"NA\"), p}' %s -",
                 cfg$cc, cfg$pc, cfg$ic, cfg$ec, cfg$oc, cfg$c4, cfg$c5, cfg$fc, rf)
  system(sprintf("%s %s 2>/dev/null | %s > %s", cfg$rd, shQuote(file.path(G, cfg$f)), awk, shQuote(out))); out }
readc <- function(path) { nm <- c("tag", "id", "ea", "oa", "c4", "c5", "eaf", "pos")
  if (file.size(path) == 0) { d <- as.data.frame(setNames(replicate(8, character(0), simplify = FALSE), nm)); d$pos <- numeric(0); return(d) }
  d <- read.table(path, header = FALSE, stringsAsFactors = FALSE, colClasses = "character", fill = TRUE, quote = "", comment.char = "")
  names(d) <- nm; d$pos <- as.numeric(d$pos); d }
prep <- function(d, cfg) {
  suppressWarnings({ d$c4 <- as.numeric(d$c4); d$c5 <- as.numeric(d$c5); d$eaf <- as.numeric(d$eaf) })
  qn <- function(p) qnorm(1 - pmin(pmax(p, 1e-300), 0.9999999) / 2)
  if (cfg$fl == "beta") { b <- d$c4; se <- d$c5 } else if (cfg$fl == "or") { b <- log(d$c4); se <- d$c5 } else { b <- d$c4; se <- abs(d$c4) / qn(d$c5) }
  d$b <- as.numeric(b); d$vb <- as.numeric(se^2); d$maf <- as.numeric(pmin(d$eaf, 1 - d$eaf))
  d[is.finite(d$b) & is.finite(d$vb) & d$vb > 0, ] }
clean <- function(d, cfg) { d <- d[grepl("^rs", d$id), ]; d <- d[!duplicated(d$id), ]; prep(d, cfg) }
mkD <- function(d, cfg) { x <- list(beta = d$b, varbeta = d$vb, snp = d$id, MAF = d$maf)
  if (cfg$ty == "quant") { x$type <- "quant"; x$N <- cfg$N } else { x$type <- "cc"; x$s <- cfg$s }; x }

## ---- 1. one pass per file: +-1 Mb around each lead variant ----
for (tr in unique(c(PAIRS$t1, PAIRS$t2))) {
  rows <- PAIRS[PAIRS$t1 == tr | PAIRS$t2 == tr, ]; b38 <- F[[tr]]$b == 38
  extract(tr, data.frame(tag = rows$id, chr = rows$chr, lo = (if (b38) rows$p38 else rows$p19) - WMAX, hi = (if (b38) rows$p38 else rows$p19) + WMAX)) }

getpair <- function(P, R, mafilter = TRUE) {   # both traits within window R, matched by rsID (OSA: hg38 position + alleles)
  # mafilter: the trait-MAF > 0.01 filter of the coloc.abf screen; the SuSiE step filters on 1000 Genomes MAF instead
  c1 <- F[[P$t1]]; c2 <- F[[P$t2]]
  w <- function(d, cfg) { ctr <- if (cfg$b == 38) P$p38 else P$p19; d[d$pos >= ctr - R & d$pos <= ctr + R, ] }
  a <- clean(w(subset(readc(file.path(CACHE, paste0(P$t1, ".tsv"))), tag == P$id), c1), c1)
  b <- w(subset(readc(file.path(CACHE, paste0(P$t2, ".tsv"))), tag == P$id), c2)
  if (P$t2 == "osa") { k <- function(x) paste(x$pos, pmin(toupper(x$ea), toupper(x$oa)), pmax(toupper(x$ea), toupper(x$oa)))
    a$pk <- k(a); b$pk <- k(b); b <- b[!duplicated(b$pk), ]; b$id <- a$id[match(b$pk, a$pk)]; b <- b[!is.na(b$id), ] }
  b <- clean(b, c2)
  cm <- intersect(a$id, b$id); A <- a[match(cm, a$id), ]; B <- b[match(cm, b$id), ]
  if (mafilter) { maf <- B$maf; maf[!is.finite(maf)] <- A$maf[!is.finite(maf)]; ok <- is.finite(maf) & maf > 0.01
    A <- A[ok, ]; B <- B[ok, ]; A$maf <- maf[ok]; B$maf <- maf[ok] }
  list(A = A, B = B, c1 = c1, c2 = c2) }

## ---- 2. window x prior ----
W <- list()
for (i in seq_len(nrow(PAIRS))) { P <- PAIRS[i, ]
  for (R in WINS) { X <- getpair(P, R); n <- nrow(X$A); if (n < 20) next
    row <- data.frame(id = P$id, window_kb = R / 1000, nSNP = n)
    for (pn in names(PRI)) { r <- suppressWarnings(coloc.abf(mkD(X$A, X$c1), mkD(X$B, X$c2), p12 = PRI[[pn]]))
      row[[paste0("PP4_", pn)]] <- round(as.numeric(r$summary["PP.H4.abf"]), 3); row[[paste0("PP3_", pn)]] <- round(as.numeric(r$summary["PP.H3.abf"]), 3)
      if (pn == "default") row$lead <- r$results$snp[which.max(r$results$SNP.PP.H4)] }
    W[[length(W) + 1]] <- row } }
W <- do.call(rbind, W)
write.table(W, file.path(RES, "coloc_window_prior_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n===== window x prior (PP.H4) =====\n"); print(W[, c("id", "window_kb", "nSNP", "PP4_strict", "PP4_default", "PP4_lenient", "PP3_default", "lead")], row.names = FALSE)

## ---- 3. coloc-SuSiE (+-250 kb) ----
ldmat <- function(snps, tag, chr, lo, hi) {
  sf <- file.path(WD, paste0(tag, "_snps.txt")); writeLines(snps, sf); sb <- file.path(WD, paste0(tag, "_sub")); lg <- file.path(WD, "plink.log")
  system2(PLINK, c("--bfile", file.path(BDIR, sprintf("1000G.EUR.QC.%d", chr)), "--extract", sf, "--chr", chr, "--from-bp", lo, "--to-bp", hi,
                   "--keep-allele-order", "--make-bed", "--out", sb), stdout = lg, stderr = lg)
  if (!file.exists(paste0(sb, ".bim"))) return(NULL)
  bim <- read.table(paste0(sb, ".bim"), stringsAsFactors = FALSE, col.names = c("chr", "rsid", "cm", "pos", "A1", "A2"))
  system2(PLINK, c("--bfile", sb, "--keep-allele-order", "--r", "square", "bin", "--out", file.path(WD, tag)), stdout = lg, stderr = lg)
  system2(PLINK, c("--bfile", sb, "--keep-allele-order", "--freq", "--out", file.path(WD, tag)), stdout = lg, stderr = lg)
  n <- nrow(bim); M <- matrix(readBin(file.path(WD, paste0(tag, ".ld.bin")), "double", n = n * n, size = 8), n, n)
  M[!is.finite(M)] <- 0; diag(M) <- 1; dimnames(M) <- list(bim$rsid, bim$rsid)
  fq <- read.table(file.path(WD, paste0(tag, ".frq")), header = TRUE, stringsAsFactors = FALSE)
  rm <- fq$MAF[match(bim$rsid, fq$SNP)]; list(LD = M, bim = bim, refmaf = pmin(rm, 1 - rm)) }
align <- function(d, bim) { d <- d[match(bim$rsid, d$id), ]
  ea <- toupper(d$ea); oa <- toupper(d$oa); a1 <- toupper(bim$A1); a2 <- toupper(bim$A2)
  same <- ea == a1 & oa == a2; flip <- ea == a2 & oa == a1
  keep <- (same | flip) & is.finite(d$b) & is.finite(d$vb) & d$vb > 0; keep[is.na(keep)] <- FALSE
  b <- d$b; b[flip & !is.na(flip)] <- -b[flip & !is.na(flip)]; list(beta = b, vb = d$vb, keep = keep) }
csdesc <- function(f) { if (is.null(f) || is.null(f$sets$cs)) return("0")
  paste0(length(f$sets$cs), ": ", paste(sapply(f$sets$cs, function(s) sprintf("%s(%d)", names(f$pip)[s][which.max(f$pip[s])], length(s))), collapse = "; ")) }

S <- list()
for (i in seq_len(nrow(PAIRS))) { P <- PAIRS[i, ]; cat(sprintf("\n--- SuSiE %s ---\n", P$id))
  X <- getpair(P, 2.5e5, mafilter = FALSE); L <- ldmat(X$A$id, P$id, P$chr, P$p19 - 2.5e5, P$p19 + 2.5e5)
  if (is.null(L)) { cat("  LD matrix could not be built\n"); next }
  fa <- align(X$A, L$bim); fb <- align(X$B, L$bim)
  ok <- fa$keep & fb$keep & is.finite(L$refmaf) & L$refmaf > 0.01; ok[is.na(ok)] <- FALSE
  cat(sprintf("  shared %d | in 1000G %d | aligned with MAF > 0.01: %d\n", nrow(X$A), nrow(L$bim), sum(ok)))
  mk <- function(f, cfg) { D <- list(beta = f$beta[ok], varbeta = f$vb[ok], snp = L$bim$rsid[ok], MAF = L$refmaf[ok], LD = L$LD[ok, ok, drop = FALSE], N = cfg$N, position = L$bim$pos[ok])
    if (cfg$ty == "quant") D$type <- "quant" else { D$type <- "cc"; D$s <- cfg$s }; D }
  FA <- tryCatch(runsusie(mk(fa, X$c1), n = X$c1$N, repeat_until_convergence = FALSE, max_iter = 300), error = function(e) NULL)
  FB <- tryCatch(runsusie(mk(fb, X$c2), n = X$c2$N, repeat_until_convergence = FALSE, max_iter = 300), error = function(e) NULL)
  row <- data.frame(id = P$id, nSNP = sum(ok), cs_t1 = csdesc(FA), cs_t2 = csdesc(FB), PP4_max = NA, hit1 = NA, hit2 = NA, PP3_of_best = NA)
  if (!is.null(FA) && !is.null(FB) && !is.null(FA$sets$cs) && !is.null(FB$sets$cs)) {
    r <- tryCatch(coloc.susie(FA, FB), error = function(e) NULL)
    if (!is.null(r) && !is.null(r$summary) && nrow(r$summary) > 0) { j <- which.max(r$summary$PP.H4.abf)
      row$PP4_max <- round(r$summary$PP.H4.abf[j], 3); row$PP3_of_best <- round(r$summary$PP.H3.abf[j], 3)
      row$hit1 <- r$summary$hit1[j]; row$hit2 <- r$summary$hit2[j]
      write.csv(r$summary, file.path(RES, paste0("susie_", P$id, ".csv")), row.names = FALSE) } }
  S[[length(S) + 1]] <- row; print(row, row.names = FALSE) }
S <- do.call(rbind, S)
write.table(S, file.path(RES, "coloc_susie_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n===== coloc-SuSiE summary =====\n"); print(S, row.names = FALSE)
cat("=== COLOC_SENSITIVITY_DONE ===\n")
