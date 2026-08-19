#!/usr/bin/env Rscript
# Ayni blok, ayni sumstats, ayni overlap — sadece LD referansi farkli.
suppressMessages(library(LAVA))
L    <- "/path/to/analysis/results/lava"
INFO <- file.path(L, "info.nalls_akcimen.txt")
OVL  <- file.path(L, "ovl.nalls_akcimen.txt")
LOCF <- file.path(L, "tox3.locfile")
G1K  <- "/path/to/analysis/data/reference/g1000_eur"
UKB  <- "/path/to/ukb_ld/lava-ukb-v1.1_chr16"
GATE <- 0.05 / 2495

loci <- read.loci(LOCF)
res <- list()
for (nm in c("1000G EUR", "UK Biobank")) {
  ref <- if (nm == "1000G EUR") G1K else UKB
  cat("\n########", nm, "########\n")
  inp <- tryCatch(process.input(INFO, OVL, ref, c("pd_nalls", "rls_akcimen")),
                  error = function(e) { cat("input err:", conditionMessage(e), "\n"); NULL })
  if (is.null(inp)) next
  loc <- tryCatch(process.locus(loci[1, ], inp), error = function(e) NULL)
  if (is.null(loc)) { cat("blok islenemedi\n"); next }
  u <- run.univ(loc); b <- run.bivar(loc)
  cat("SNP:", loc$n.snps, " | bilesen (K):", loc$K, "\n")
  print(u)
  cat("\n")
  print(b)
  res[[nm]] <- data.frame(
    reference = nm, n.snps = loc$n.snps, K = loc$K,
    h2_pd  = signif(u$h2.obs[u$phen == "pd_nalls"], 3),
    p_univ_pd = signif(u$p[u$phen == "pd_nalls"], 3),
    gate_pd = u$p[u$phen == "pd_nalls"] < GATE,
    h2_rls = signif(u$h2.obs[u$phen == "rls_akcimen"], 3),
    p_univ_rls = signif(u$p[u$phen == "rls_akcimen"], 3),
    gate_rls = u$p[u$phen == "rls_akcimen"] < GATE,
    rho = round(b$rho[1], 4), lo = round(b$rho.lower[1], 4),
    hi = round(b$rho.upper[1], 4), p_bivar = signif(b$p[1], 3))
}
cat("\n\n======== YAN YANA ========\n")
out <- do.call(rbind, res)
print(out, row.names = FALSE)
cat("\nunivariate kapi: p <", signif(GATE, 3), "(her iki traitte de gecmeli)\n")

# UKB referansinin bu bolgedeki SNP kapsami
info <- try(read.table("/path/to/ukb_ld/lava-ukb-v1.1_chr16.info", header = TRUE,
                       nrows = 5), silent = TRUE)
if (!inherits(info, "try-error")) {
  cat("\n.info kolonlari:", paste(names(info), collapse = ", "), "\n")
}
write.csv(out, "/path/to/analysis/results/lava/tox3_ld_reference_comparison.csv",
          row.names = FALSE)
cat("\nyazildi: results/lava/tox3_ld_reference_comparison.csv\n")
