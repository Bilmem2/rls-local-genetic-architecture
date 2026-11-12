# Results

Machine-readable primary results behind the manuscript.

- **`global_genetic_correlation_ldsc.csv`** - LDSC genome-wide genetic correlation of RLS (Akcimen)
  and Parkinson's disease (Nalls) with each panel trait (`rg`, standard error, z, p, and the
  cross-trait LDSC `gcov_intercept` that accounts for sample overlap). Backs manuscript Supplementary Table S2:
  RLS correlates broadly with sleep and internalizing traits, PD only with essential tremor.
- **`colocalization_hits.csv`** - the colocalizing loci (`coloc.abf` PP.H4 > 0.90) plus the single
  suggestive locus (0.80-0.90), with the implicated gene and lead variant (manuscript Table 1).
- **`tox3_crossancestry_local_rg.csv`** - TOX3/CASC16 (chr16q12) local genetic correlation between
  RLS and PD across PD GWAS and ancestries (EUR/AFR/AMR); all estimates negative (manuscript Figure 5).
- **`xanc_coloc_afr.csv`** - African-ancestry local genetic correlation re-estimated at each
  European colocalizing locus for which African-ancestry data exist for both traits. The `TOX3` row
  reproduces the primary cross-ancestry estimate above (rho = -0.734) and so acts as a positive
  control on this screen; the sleep loci are wide and non-significant, which reflects the size of
  the African-ancestry RLS arm (2,176 cases) rather than a failure to replicate, and `SEMA6D` is
  not estimable because RLS shows no detectable African-ancestry local heritability there
  (manuscript Discussion).
- **`coloc_screen_full_supplementary.csv`** - every local genetic correlation reaching the LAVA
  screening threshold (p < 1e-3) carried to colocalization, with its coloc result. Also carries
  the multiple-testing columns: `p_bonferroni` (over all 2,821 bivariate tests admitted by the
  both-trait local-heritability gate), `q_BH` (Benjamini-Hochberg over the same tests) and
  `survives_bonferroni`. All 74 rows are FDR-significant; 14 also survive Bonferroni
  (manuscript Supplementary Table S3).
- **`known_vs_novel_blocks.csv`** - for each of the 74 screened local correlations, the smallest P and its
  variant in each trait's GWAS within the block, and a class: `both`, `trait 1 only`, `trait 2 only`
  or `neither` carry a genome-wide-significant variant (P < 5e-8). Of the 71 assessable rows, 10 are
  `both`, 14 `trait 1 only`, 8 `trait 2 only` and 39 `neither`; all seven colocalizing pairs carry
  a genome-wide-significant variant in at least one trait. The three narcolepsy rows cannot be
  assessed: the public Ollila release removes the MHC region and every variant with P < 5e-8.
  Block membership is by rsID against the 1000 Genomes EUR reference, so hg38 and hg19 files are
  treated alike (GWS column of manuscript Supplementary Table S3).
- **`coloc_prior_sensitivity.tsv`** - PP.H4 for each prioritized locus recomputed across
  colocalization priors (p12 from 1e-6 to 5e-5, against the coloc default of 1e-5). TOX3, MEIS1
  and SKOR2 hold under the strictest prior; PTPRD falls to 0.51 and the suggestive GBA signal to
  0.31 (manuscript Supplementary Table S4, panel A).
- **`coloc_window_sensitivity.tsv`** - PP.H4 recomputed over windows of +-100 kb to +-1 Mb
  centred on each lead variant. Flat everywhere except SEMA6D, which is lost at +-1 Mb where a
  distinct neighbouring association enters (Supplementary Table S4, panel B).
- **`susie_meis1_rls_*.csv`** - coloc-SuSiE credible-set pairings at MEIS1 for RLS against
  insomnia, chronotype and short sleep. RLS resolves into 8-10 independent credible sets, yet
  each sleep trait colocalizes with exactly one of them (Supplementary Table S4, panel C).
- **`susie_tox3_rls_pd_nalls.csv`, `susie_sema6d_rls_mdd.csv`** - coloc-SuSiE credible-set pairings at
  TOX3 (RLS x PD) and SEMA6D (RLS x MDD). At TOX3 RLS resolves into two credible sets and only the one
  led by rs3104769 colocalizes with PD (PP.H4 0.87); at SEMA6D each trait has one set (PP.H4 0.86).
  Both are suggestive (0.80-0.90) by the manuscript's thresholds (Supplementary Table S4, panel C).
- **`coloc_susie_summary.tsv`** - every credible-set pairing at the five colocalizing loci, with
  `coloc.abf` PP.H4 on the same SNP set. At PTPRD and SKOR2 SuSiE finds no credible set for the
  sleep trait, so coloc-SuSiE cannot be applied (`coloc.abf` on the same variants: 0.91 and 0.98).
- **`lava_pd_rls_afr.csv`, `lava_pd_rls_amr.csv`** - genome-wide gated LAVA scans of RLS x PD in
  African (272 blocks) and Admixed-American (196 blocks) ancestry. TOX3 is block 2134: in the African
  scan P=2.7e-5, Bonferroni-corrected 7.4e-3 (`p_bonf`). The block did not pass the local-heritability
  gate in the Admixed-American scan, so the AMR estimate in Figure 5 is an ungated point estimate.
- **`functional_annotation.tsv`** - manuscript Table 2: QTL, credible-set, L2G and enhancer evidence
  for each colocalizing lead variant, taken by hand from the Open Targets Platform (release 26.06,
  July 2026) and GTEx.

Numbers are taken directly from the analysis outputs (LDSC logs, the coloc screen, and the gated
LAVA estimates); see `scripts/` for how each was produced.
