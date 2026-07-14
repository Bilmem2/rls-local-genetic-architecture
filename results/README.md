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
- **`xanc_coloc_afr.csv`** - African-ancestry local genetic correlation re-estimated at each
  European colocalizing locus for which African-ancestry data exist for both traits. The `TOX3` row
  reproduces the primary cross-ancestry estimate above (rho = -0.734) and so acts as a positive
  control on this screen; the sleep loci are wide and non-significant, which reflects the size of
  the African-ancestry RLS arm (2,176 cases) rather than a failure to replicate, and `SEMA6D` is
  not estimable because RLS shows no detectable African-ancestry local heritability there
  (manuscript Discussion).
- **`coloc_screen_full_supplementary.csv`** - every significant local genetic correlation (LAVA
  p < 1e-3) carried to colocalization, with its coloc result (Supplementary Table S1).

Numbers are taken directly from the analysis outputs (LDSC logs, the coloc screen, and the gated
LAVA estimates); see `scripts/` for how each was produced.
