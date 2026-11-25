#!/usr/bin/env Rscript
# CCDC148/PKP4 (chr2q24): RLS (Akcimen 2026, European arm) x obstructive sleep apnea from FinnGen R12
# (G6_SLEEPAPNO; 56,885 cases / 441,137 controls; hg38). FinnGen is not part of the European RLS
# meta-analysis (All of Us, UK Biobank, CLSA, CARTaGENE, MVP), so the two share no participants; FinnGen
# is, however, one of the cohorts of the discovery OSA meta-analysis.
# Same analyses as for the discovery OSA data: coloc.abf over the LAVA block under three priors and over
# four windows around the lead variant (rs77571593), and coloc-SuSiE over +-250 kb with 1000 Genomes EUR LD.
# Output: results/finngen_osa_coloc.tsv (Supplementary Table S4, panel E).
suppressMessages({ library(coloc); library(susieR) }); options(scipen = 999)
ROOT <- "/path/to/analysis"; G <- file.path(ROOT, "data/gwas")
REPO <- Sys.getenv("REPO", "."); OUT <- file.path(REPO, "results/finngen_osa_coloc.tsv")
WD <- file.path(ROOT, "results/susie_ld"); dir.create(WD, showWarnings = FALSE, recursive = TRUE)
PLINK <- Sys.getenv("PLINK", "plink"); BFILE <- file.path(ROOT, "data/ldsc_ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.2")
RLSF <- file.path(G, "rls/akcimen2026_EUR_RLS.txt"); FIN <- file.path(G, "sleep/finngen_osa_r12/finngen_R12_G6_SLEEPAPNO.gz")
S3 <- read.csv(file.path(REPO, "results/coloc_screen_full_supplementary.csv"), stringsAsFactors = FALSE)
r <- S3[S3$Trait_1 == "RLS" & S3$Trait_2 == "OSA" & S3$Region_hg19 == "chr2:159040980-160200379", ]; stopifnot(nrow(r) == 1)
B38 <- as.numeric(strsplit(r$region_hg38, "[:-]")[[1]]); LEAD38 <- 158443890; LEAD19 <- 159300402
cat(sprintf("LAVA block hg38: chr%d:%d-%d | lead rs77571593 hg38 %d\n", B38[1], B38[2], B38[3], LEAD38))

hdr <- strsplit(system(sprintf("zcat %s | awk 'NR==1'", shQuote(FIN)), intern = TRUE), "\t")[[1]]
stopifnot(identical(hdr[c(1:4, 7, 9:11)], c("#chrom", "pos", "ref", "alt", "pval", "beta", "sebeta", "af_alt")))

lo <- min(B38[2], LEAD38 - 1e6); hi <- max(B38[3], LEAD38 + 1e6)
rls <- read.table(pipe(sprintf("awk '$2==2&&$3>=%.0f&&$3<=%.0f&&$10~/^rs/{print $10,$3,$4,$5,$6,$7,$8}' %s", lo, hi, shQuote(RLSF))),
                  col.names = c("id", "pos", "ea", "oa", "b", "se", "p"), stringsAsFactors = FALSE)
rls <- rls[!duplicated(rls$id), ]
fin <- read.table(pipe(sprintf("zcat %s | awk 'NR>1&&($1==\"2\"||$1==\"chr2\")&&$2>=%.0f&&$2<=%.0f{print $2,$3,$4,$7,$9,$10,$11}'", shQuote(FIN), lo, hi)),
                  col.names = c("pos", "ref", "alt", "p", "b", "se", "af"), stringsAsFactors = FALSE)
key <- function(p, a, b) paste(p, pmin(toupper(a), toupper(b)), pmax(toupper(a), toupper(b)))
rls$k <- key(rls$pos, rls$ea, rls$oa); fin$k <- key(fin$pos, fin$ref, fin$alt); fin <- fin[!duplicated(fin$k), ]
m <- merge(rls, fin, by = "k", suffixes = c(".r", ".f"))
m$bf <- ifelse(toupper(m$alt) == toupper(m$ea), m$b.f, -m$b.f)     # FinnGen effect aligned to the RLS effect allele
m$maf <- pmin(m$af, 1 - m$af); m <- m[is.finite(m$bf) & is.finite(m$se.f) & m$se.f > 0 & is.finite(m$se.r) & m$se.r > 0 & m$maf > 0.01, ]
cat(sprintf("RLS %d | FinnGen %d | matched on position and alleles (MAF > 0.01): %d\n", nrow(rls), nrow(fin), nrow(m)))
L <- m[m$id == "rs77571593", ]
cat(sprintf("rs77571593: RLS beta %.4f P %.2g | FinnGen OSA beta (same allele) %.4f P %.2g\n", L$b.r, L$p.r, L$bf, L$p.f))
inb <- m$pos.r >= B38[2] & m$pos.r <= B38[3]
cat(sprintf("Smallest FinnGen OSA P in the block: %.2g (%s)\n", min(m$p.f[inb]), m$id[inb][which.min(m$p.f[inb])]))

D1 <- function(x) list(beta = x$b.r, varbeta = x$se.r^2, snp = x$id, MAF = x$maf, type = "cc", s = 0.0547)
D2 <- function(x) list(beta = x$bf, varbeta = x$se.f^2, snp = x$id, MAF = x$maf, type = "cc", s = 56885 / 498022)
res <- list()
cc <- function(x, lab) { row <- data.frame(analysis = lab, nSNP = nrow(x))
  for (pn in c(strict = 1e-6, default = 1e-5, lenient = 5e-5)) { rr <- suppressWarnings(coloc.abf(D1(x), D2(x), p12 = pn))
    nm <- names(which(c(strict = 1e-6, default = 1e-5, lenient = 5e-5) == pn))
    row[[paste0("PP4_", nm)]] <- round(as.numeric(rr$summary["PP.H4.abf"]), 3); row[[paste0("PP3_", nm)]] <- round(as.numeric(rr$summary["PP.H3.abf"]), 3)
    if (pn == 1e-5) row$lead <- rr$results$snp[which.max(rr$results$SNP.PP.H4)] }
  row }
res[[1]] <- cc(m[inb, ], "LAVA block")
for (R in c(1e5, 2.5e5, 5e5, 1e6)) res[[length(res) + 1]] <- cc(m[abs(m$pos.r - LEAD38) <= R, ], sprintf("lead +-%d kb", R / 1000))
W <- do.call(rbind, res); cat("\n===== RLS x FinnGen OSA: coloc.abf =====\n"); print(W, row.names = FALSE)

## coloc-SuSiE (+-250 kb), settings as in coloc_sensitivity.R
x <- m[abs(m$pos.r - LEAD38) <= 2.5e5, ]; sf <- file.path(WD, "fin_snps.txt"); writeLines(x$id, sf); sb <- file.path(WD, "fin_sub"); lg <- file.path(WD, "plink.log")
system2(PLINK, c("--bfile", BFILE, "--extract", sf, "--chr", 2, "--from-bp", LEAD19 - 2.5e5, "--to-bp", LEAD19 + 2.5e5, "--keep-allele-order", "--make-bed", "--out", sb), stdout = lg, stderr = lg)
bim <- read.table(paste0(sb, ".bim"), col.names = c("chr", "rsid", "cm", "pos", "A1", "A2"), stringsAsFactors = FALSE)
system2(PLINK, c("--bfile", sb, "--keep-allele-order", "--r", "square", "bin", "--out", file.path(WD, "fin")), stdout = lg, stderr = lg)
system2(PLINK, c("--bfile", sb, "--keep-allele-order", "--freq", "--out", file.path(WD, "fin")), stdout = lg, stderr = lg)
n <- nrow(bim); LD <- matrix(readBin(file.path(WD, "fin.ld.bin"), "double", n = n * n, size = 8), n, n); LD[!is.finite(LD)] <- 0; diag(LD) <- 1
dimnames(LD) <- list(bim$rsid, bim$rsid); fq <- read.table(file.path(WD, "fin.frq"), header = TRUE); rm <- fq$MAF[match(bim$rsid, fq$SNP)]; rm <- pmin(rm, 1 - rm)
y <- x[match(bim$rsid, x$id), ]; same <- toupper(y$ea) == toupper(bim$A1) & toupper(y$oa) == toupper(bim$A2); flip <- toupper(y$ea) == toupper(bim$A2) & toupper(y$oa) == toupper(bim$A1)
ok <- (same | flip) & is.finite(rm) & rm > 0.01; ok[is.na(ok)] <- FALSE; sg <- ifelse(flip, -1, 1)
cat(sprintf("\nSuSiE: shared %d | in 1000G %d | aligned with MAF > 0.01: %d\n", nrow(x), n, sum(ok)))
mk <- function(b, se, N, s) list(beta = (sg * b)[ok], varbeta = (se^2)[ok], snp = bim$rsid[ok], MAF = rm[ok], LD = LD[ok, ok], N = N, position = bim$pos[ok], type = "cc", s = s)
FA <- tryCatch(runsusie(mk(y$b.r, y$se.r, 676175, 0.0547), n = 676175, repeat_until_convergence = FALSE, max_iter = 300), error = function(e) NULL)
FB <- tryCatch(runsusie(mk(y$bf, y$se.f, 498022, 56885 / 498022), n = 498022, repeat_until_convergence = FALSE, max_iter = 300), error = function(e) NULL)
cs <- function(f) if (is.null(f) || is.null(f$sets$cs)) "0" else paste(length(f$sets$cs), paste(sapply(f$sets$cs, function(s) sprintf("%s(%d)", names(f$pip)[s][which.max(f$pip[s])], length(s))), collapse = "; "))
cat("RLS credible sets:", cs(FA), "\nFinnGen OSA credible sets:", cs(FB), "\n")
if (!is.null(FA$sets$cs) && !is.null(FB$sets$cs)) { cr <- coloc.susie(FA, FB); print(cr$summary[, c("hit1", "hit2", "PP.H3.abf", "PP.H4.abf")], row.names = FALSE) } else cat("coloc-SuSiE not applicable: no credible set for one trait\n")
write.table(W, OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== FINNGEN_OSA_DONE ===\n")
