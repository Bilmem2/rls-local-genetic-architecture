# Results

Machine-readable primary results behind the manuscript.

- **`global_genetic_correlation_ldsc.csv`** - LDSC genome-wide genetic correlation of RLS (Akcimen)
  and Parkinson's disease (Nalls) with each panel trait (`rg`, standard error, z, p, and the
  cross-trait LDSC `gcov_intercept` that accounts for sample overlap). Backs manuscript Supplementary Table S2:
  RLS correlates broadly with sleep and internalizing traits, PD only with essential tremor.
- **`lava_test_family.tsv`** - every bivariate LAVA test of the European panel admitted by the
  both-trait local-heritability gate, each test once (2,739 tests), with `rho`, `p`, `p_bonf`
  (Bonferroni), `q_BH` (Benjamini-Hochberg) and the LAVA output file it was taken from. 74 tests reach
  the screening threshold (P < 1e-3), all with q_BH < 0.05; 14 survive Bonferroni (P < 1.8e-5).
- **`colocalization_hits.csv`** - the seven colocalizing pairs (`coloc.abf` PP.H4 > 0.90) plus the
  suggestive one (0.80-0.90), with the implicated gene, the lead variant and the evidence label of
  manuscript Table 1: `robust` pairs keep PP.H4 >= 0.80 under every prior and window tested;
  `prior-sensitive` pairs fall below 0.80 under p12 = 1e-6, `window-sensitive` pairs under the wider windows.
- **`coloc_screen_full_supplementary.csv`** - every local genetic correlation reaching the LAVA
  screening threshold, with its `coloc.abf` result over the LAVA block (PP.H3 and PP.H4 under the
  default prior, PP.H4 under p12 = 1e-6 and 5e-5), the lead variant, the verdict (`COLOC` PP.H4 > 0.90,
  `suggestive` 0.80-0.90), the Bonferroni and BH columns of `lava_test_family.tsv` and the GWS class
  from `known_vs_novel_blocks.csv`. Seven rows colocalize and one is suggestive (manuscript
  Supplementary Table S3).
- **`known_vs_novel_blocks.csv`** - for each of the 74 screened local correlations, the smallest P and its
  variant in each trait's GWAS within the block, and a class: `both`, `trait 1 only`, `trait 2 only`
  or `neither` carry a genome-wide-significant variant (P < 5e-8). Of the 71 assessable rows, 10 are
  `both`, 14 `trait 1 only`, 8 `trait 2 only` and 39 `neither`; all eight colocalizing or suggestive
  pairs carry a genome-wide-significant variant in at least one trait. The three narcolepsy rows cannot
  be assessed: the public Ollila release removes the MHC region and every variant with P < 5e-8.
  Block membership is by rsID against the 1000 Genomes EUR reference, except for obstructive sleep
  apnea (no rsIDs), where it is by position within the block's hg38 coordinates (GWS column of
  Supplementary Table S3).
- **`coloc_window_prior_sensitivity.tsv`** - PP.H4 and PP.H3 for each of the eight pairs over windows
  of +-100 kb to +-1 Mb centred on the lead variant, each under p12 = 1e-6, 1e-5 and 5e-5. *MEIS1*
  (insomnia, short sleep, chronotype) and *SKOR2* stay at or above 0.80 throughout; *TOX3* (0.56) and
  *CCDC148/PKP4* (0.48) fall below 0.80 under p12 = 1e-6; *SEMA6D* is lost at +-1 Mb, where a distinct
  neighbouring association enters (Supplementary Table S4, panels A and B).
- **`susie_*.csv`, `coloc_susie_summary.tsv`** - coloc-SuSiE over +-250 kb around each lead variant
  (LD from the 1000 Genomes Phase 3 European samples distributed with LDSC, n = 489). At *MEIS1* RLS
  resolves into 9-10 credible sets, and each sleep trait colocalizes with the set led by rs113851554
  (PP.H4 0.98-1.00). At *TOX3* RLS resolves into two credible sets and only the one led by rs3104769
  colocalizes with PD (PP.H4 0.87); *SEMA6D* (0.86) and *CCDC148/PKP4* (0.80) are suggestive by the
  manuscript's thresholds. At *SKOR2* and *GBA1* SuSiE finds no credible set for the second trait, so
  coloc-SuSiE cannot be applied (Supplementary Table S4, panel C).
- **`tox3_ld_reference_comparison.csv`** - the TOX3 RLS x PD local correlation (Nalls x Akcimen) with the
  1000 Genomes EUR and a UK Biobank LD reference, with the range of LAVA's P over eight seeds. With
  UK Biobank LD the estimate stays negative but PD local heritability no longer passes the gate
  (Supplementary Table S4, panel D).
- **`finngen_osa_coloc.tsv`** - the *CCDC148/PKP4* RLS x obstructive sleep apnea signal tested with
  FinnGen R12 (G6_SLEEPAPNO), which shares no participants with the RLS data but is one of the cohorts
  of the discovery OSA meta-analysis: PP.H4 0.23 over the LAVA block (0.03 to 0.60 across priors) and
  0.23 to 0.28 across windows; FinnGen has no credible set there, so coloc-SuSiE cannot be applied
  (Supplementary Table S4, panel E).
- **`tox3_crossancestry_local_rg.csv`** - TOX3/CASC16 (chr16q12) local genetic correlation between
  RLS and PD across PD GWAS and ancestries (EUR/AFR/AMR); all 11 estimates negative (manuscript Figure 5).
  Rows with a confidence interval are gated LAVA estimates; the others are ungated point estimates.
  The GP2 rows use the effective sample size and the LDSC-based overlap term of `xanc_lava_tox3.sh`.
- **`xanc_coloc_afr.csv`** - African-ancestry local genetic correlation re-estimated at each
  European colocalizing locus for which African-ancestry data exist for both traits. The `TOX3` row
  reproduces the primary cross-ancestry estimate (rho = -0.700) and so acts as a positive control on
  this check; the *MEIS1* short-sleep estimate is wide and non-significant, which reflects the size of
  the African-ancestry RLS arm (2,176 cases), and `SEMA6D` is not estimable because RLS shows no
  detectable African-ancestry local heritability there (manuscript Discussion).
- **`lava_pd_rls_afr.csv`, `lava_pd_rls_amr.csv`** - genome-wide gated LAVA scans of RLS x PD in
  African (314 blocks) and Admixed-American (196 blocks) ancestry. TOX3 is block 2134: in the African
  scan P=2.9e-5, Bonferroni-corrected 9.0e-3 (`p_bonf`). The block did not pass the local-heritability
  gate in the Admixed-American scan, so the AMR estimate in Figure 5 is an ungated point estimate.
- **`functional_annotation.tsv`** - manuscript Table 2: credible-set, QTL, L2G and enhancer evidence
  for each colocalizing lead variant, from the Open Targets Platform (release 26.09, November 2026).

Numbers are taken directly from the analysis outputs (LDSC logs, the coloc screen, and the gated
LAVA estimates); see `scripts/` for how each was produced.
