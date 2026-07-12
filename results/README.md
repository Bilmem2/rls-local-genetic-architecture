# Results

Machine-readable primary results behind the manuscript.

- **`global_genetic_correlation_ldsc.csv`** - LDSC genome-wide genetic correlation of RLS (Akcimen)
  and Parkinson's disease (Nalls) with each panel trait (`rg`, standard error, z, p, and the
  cross-trait LDSC `gcov_intercept` that accounts for sample overlap). Backs the "RLS is a global
  hub / PD is globally isolated" result (manuscript Table 2).
- **`colocalization_hits.csv`** - the colocalizing loci (`coloc.abf` PP.H4 > 0.90) plus the single
  suggestive locus (0.80-0.90), with the implicated gene and lead variant (manuscript Table 3).
- **`tox3_crossancestry_local_rg.csv`** - TOX3/CASC16 (chr16q12) local genetic correlation between
  RLS and PD across PD GWAS and ancestries (EUR/AFR/AMR); all estimates negative (manuscript Figure 5).
- **`coloc_screen_full_supplementary.csv`** - every significant local genetic correlation (LAVA
  p < 1e-3) carried to colocalization, with its coloc result (Supplementary Table S1).

Numbers are taken directly from the analysis outputs (LDSC logs, the coloc screen, and the gated
LAVA estimates); see `scripts/` for how each was produced.
