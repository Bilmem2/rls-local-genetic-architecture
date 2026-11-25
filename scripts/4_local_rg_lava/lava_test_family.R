#!/usr/bin/env Rscript
# Unique LAVA bivariate test family for the European panel, with Bonferroni and Benjamini-Hochberg
# adjustment (p.adjust).
#
# Inputs are the bivariate LAVA outputs of panel_lava*.R and run_narco_lava.R. Some tests were run in
# more than one output file; each test is kept once, from the most recently written file (a rule fixed
# without reference to the P values). Narcolepsy enters as the Ollila 2023 GWAS (narcolepsy_ollila);
# the reversed-order RLS x PD run (lava_pd_rls.csv) is not read.
#
# Output: results/lava_test_family.tsv - 2,739 tests; 74 reach the screening threshold (P < 1e-3), all
# with BH q < 0.05, and 14 survive Bonferroni (manuscript Results and Supplementary Table S3).
ROOT <- "/path/to/analysis"; LV <- file.path(ROOT, "results/lava")
REPO <- Sys.getenv("REPO", "."); OUT <- Sys.getenv("OUT", file.path(REPO, "results/lava_test_family.tsv"))
files <- c("lava_panel_bivariate.csv", "lava_panel_rlssig.csv", "lava_panel_rlstail.csv",
           "lava_panel_tail.csv", "lava_panel_tier3.csv", "lava_narco_ollila.csv")
panel <- c("anx","insomnia","mdd","osa","ptsd","rbd","sleepdur_eur_short","sleepdur_eur_long","bip","chronotype",
           "daytime_sleepiness","scz","et","tourette","narcolepsy_ollila","pd_nalls")
mt <- file.mtime(file.path(LV, files)); names(mt) <- files
cat("File modification times:\n"); print(sort(mt))
D <- do.call(rbind, lapply(files, function(f) {
  d <- read.csv(file.path(LV, f), stringsAsFactors = FALSE)
  d <- d[is.finite(d$p) & d$phen1 %in% c("rls_akcimen","pd_nalls") & d$phen2 %in% panel, c("phen1","phen2","chr","start","stop","rho","p")]
  d$file <- f; d$mtime <- mt[f]; d }))
D <- D[!(D$phen1 == "pd_nalls" & D$phen2 == "pd_nalls"), ]
D$key <- paste(D$phen1, D$phen2, D$chr, D$start)
cat(sprintf("\nRows read: %d | unique tests: %d\n", nrow(D), length(unique(D$key))))
dupk <- unique(D$key[duplicated(D$key)])
cat(sprintf("Tests present in more than one file: %d\n", length(dupk)))
D <- D[order(D$key, -as.numeric(D$mtime)), ]
U <- D[!duplicated(D$key), ]
stopifnot(!any(duplicated(U$key)))
cat("File kept for the repeated tests:\n"); print(table(U$file[U$key %in% dupk]))
m <- nrow(U)
U$p_bonf <- p.adjust(U$p, "bonferroni"); U$q_BH <- p.adjust(U$p, "BH")
stopifnot(all(U$q_BH <= U$p_bonf + 1e-15))   # BH can never exceed Bonferroni
cat(sprintf("\nTest family m = %d | Bonferroni threshold 0.05/m = %.3e\n", m, 0.05 / m))
cat("Tests per trait pair:\n"); print(table(paste(U$phen1, U$phen2)))
H <- U[U$p < 1e-3, ]
cat(sprintf("\nReach the screening threshold (P < 1e-3): %d | survive Bonferroni: %d | all BH q < 0.05: %s\n",
            nrow(H), sum(U$p_bonf < 0.05), all(H$q_BH < 0.05)))
write.table(U[, c("phen1","phen2","chr","start","stop","rho","p","p_bonf","q_BH","file")],
            OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("written:", OUT, "\n")
