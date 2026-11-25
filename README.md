# Local genetic architecture of Restless Legs Syndrome

Analysis code for the manuscript **"Local Genetic Architecture of Restless Legs Syndrome and
Parkinson's Disease Across Sleep, Psychiatric, and Movement Disorders."**

The pipeline is anchored on Restless Legs Syndrome (RLS) and screens
a panel of sleep, psychiatric, and movement traits, including Parkinson's disease (PD), through
four layers:

1. **Global** genetic correlation (LDSC)
2. **Local** genetic correlation (LAVA, per ~2,495 LD blocks)
3. **Colocalization** (`coloc.abf` at every block passing the local screen, plus prior/window
   sensitivity and `coloc.susie` at the colocalizing pairs)
4. **Cross-ancestry** comparison of the TOX3 RLS-PD signal (EUR / AFR / AMR)

Summary-statistic files and LD reference panels are **not** redistributed here (they are public or
collaborator-provided; see **Data**). The scripts document every input and its source.

## Running the pipeline

After setting the paths (see **Configuration**), the top-level driver runs every stage in order:

    bash run_all.sh

`run_all.sh` is the single source of truth for execution order. In brief:

1. **Setup** - build the toolchain; download LD references (1000G EUR/AFR/AMR, HapMap3).
2. **Munge** - download and harmonize every GWAS (RLS, PD, panel).
3. **Global rg** - LDSC: RLS x panel and PD x panel.
4. **Local rg** - LAVA genome-wide, then cross-ancestry LAVA at TOX3.
5. **Coloc** - systematic `coloc.abf` screen, then the robustness analyses (prior, window,
   coloc-SuSiE, FinnGen obstructive sleep apnea).
6. **Cross-ancestry** - TOX3 forest + per-SNP (beta-beta) concordance.
7. **Figures**.

Within a stage the scripts are largely independent, with two ordering rules that `run_all.sh`
already encodes: each `annotate_*` script must run before the `munge_*` scripts that consume its
rsID-annotated output, and the coloc / figure steps consume the LAVA output. These are the scripts
as run; the RLS-anchored steps require the access-controlled multi-ancestry RLS data (see
**Reproducibility**).

## Repository layout

Scripts are grouped by pipeline stage; run the stages in order (or use `run_all.sh`).

### `scripts/1_setup/` - toolchain and LD references
- `00_setup_environment.sh` - build the version-pinned LDSC / LAVA / PLINK toolchain (records the LDSC Python-3 patches)
- `01b_get_ldsc_ref.sh` - LDSC LD scores + HapMap3 weights (1000G Phase 3)
- `01c_smoketest.sh` - end-to-end LDSC validation on synthetic data
- `04_get_1000g_eur.sh` - 1000G Phase 3 EUR PLINK panel (LAVA / clumping)
- `dl_g1000_panels.sh` - 1000G AFR + AMR panels (cross-ancestry LAVA)
- `build_hm3_hg38_map.sh` - HapMap3 hg38->rsID map (to annotate hg38 PD sumstats)

### `scripts/2_munge/` - download and harmonize all GWAS
- RLS: `16_munge_akcimen.sh` (Akcimen EUR), `dl_rls_xanc.sh` + `munge_rls_xanc.sh` (multi-ancestry AFR/AMR)
- PD: `pd_01_download.sh` + `pd_02_munge.sh` (Nalls + FinnGen), `munge_pd_rizig.sh` (AFR), `munge_pd_xanc.sh` + `annotate_pd_xanc.py` (GP2 EUR / Loesch AMR), `munge_pd_largepd.sh` (LARGE-PD AMR)
- Panel: `annotate_sleep.sh` (rsID annotation), `munge_panel.sh` (psychiatric + sleep), `munge_insomnia.sh`, `munge_sleep2.sh` (OSA + sleep duration), `munge_sleep3.sh` (sleep-duration rsID fix), `munge_rg_daytime.sh` (daytime sleepiness), `munge_narcolepsy_ollila.sh` (narcolepsy, Ollila 2023 EUR), `dl_finngen_osa.sh` (FinnGen R12 obstructive sleep apnea)

### `scripts/3_global_rg/` - LDSC genetic correlation
- `panel_rg_akcimen.sh` (primary: RLS-Akcimen x panel), `panel_rg.sh` (RLS + PD x panel), `insomnia_rg.sh`, `pd_03_rg.sh` (PD x RLS)

### `scripts/4_local_rg_lava/` - LAVA local genetic correlation
- `panel_lava.R` (genome-wide bivariate local rg; `panel_lava2.R` is the parallelizable variant), `panel_lava_run.sh` / `panel_lava2_run.sh` (runners)
- `pd_05_lava.R` (the genome-wide LAVA engine driven by `LAVA_PAIRS`), `pd_xanc_run_lava.sh` (Admixed-American RLS x PD scan), `xanc_lava_tox3.sh` (African-ancestry RLS x PD scan, and the GP2 European meta-analysis against the three European RLS GWAS at TOX3 with the LDSC-based sample-overlap term), `neff_from_se.sh` (effective sample size of the GP2 meta-analysis from its standard errors)
- `run_narco_lava.R` (narcolepsy against RLS and PD, Ollila 2023)
- `lava_test_family.R` (the unique family of 2,739 bivariate tests, with Bonferroni and Benjamini-Hochberg adjustment)
- `known_vs_novel_blocks.R` and `gws_full_files.R` (for each of the 74 screened local correlations, whether either trait's GWAS carries a genome-wide-significant variant in the block; the second script re-reads short sleep, long sleep and obstructive sleep apnea from the full summary statistics)
- `tox3_ukb_ld.sh`, `tox3_compare_refs.R` and `lava_stochasticity.R` (TOX3 RLS x PD local correlation with a UK Biobank LD reference instead of 1000 Genomes, and the seed-to-seed variation of LAVA's P and CI; Supplementary Table S4, panel D)

### `scripts/5_coloc/` - colocalization (coloc.abf)
- `coloc_screen_full.R` (`coloc.abf` over the LAVA block of each of the 74 local correlations with P < 1e-3, every GWAS read in its own genome build from the full summary statistics, under three priors) and `build_table_s3.py` (assembles `results/coloc_screen_full_supplementary.csv`, **Supplementary Table S3**)
- robustness: `coloc_sensitivity.R` (PP.H4 over four windows around each lead variant under three priors, and coloc-SuSiE allowing several causal variants per trait; at *SKOR2* and *GBA1* SuSiE finds no credible set for the second trait, so no set-level comparison is possible) and `finngen_osa_coloc.R` (the *CCDC148/PKP4* RLS x obstructive sleep apnea signal tested in FinnGen R12, which shares no participants with the RLS data but is one of the cohorts of the discovery obstructive sleep apnea meta-analysis) - together **Supplementary Table S4**

### `scripts/6_cross_ancestry_tox3/` - TOX3/CASC16 cross-ancestry
- `run_xanc_coloc_afr.R` (African-ancestry local rg at every colocalizing locus; the TOX3 row reproduces the primary estimate and acts as a positive control)
- `forest_tox3.py` (local-rg forest across datasets and ancestries; Figure 5), `betabeta_extract.sh` + `betabeta_plot.py` (per-SNP effect concordance)

### `scripts/7_figures/`
- `forest_rg.R` (global-rg forest, Figure 2), `coloc_selectivity.py` (local rg against colocalization posterior, Figure 3)

Both read their numbers from the tables in `results/`, so a figure cannot drift away from the reported values.

### Functional annotation (no script)
The lead variants of the colocalizing pairs were annotated in the Open Targets Platform (release 26.09,
accessed November 2026): GWAS and molecular-QTL credible sets, locus-to-gene (L2G) scores and
enhancer-to-gene predictions. The extracted values are in `results/functional_annotation.tsv`
(manuscript Table 2).

## Results

`results/` holds the machine-readable numbers behind the manuscript; `results/README.md` describes
each file. In brief: `global_genetic_correlation_ldsc.csv` (LDSC, **Supplementary Table S2**),
`lava_test_family.tsv` (the 2,739 bivariate LAVA tests with Bonferroni and Benjamini-Hochberg
adjustment), `colocalization_hits.csv` (**Table 1**), `coloc_screen_full_supplementary.csv` (all 74
local correlations carried to colocalization, **Supplementary Table S3**),
`coloc_window_prior_sensitivity.tsv`, the coloc-SuSiE outputs `susie_*.csv` summarised in
`coloc_susie_summary.tsv`, `tox3_ld_reference_comparison.csv` and `finngen_osa_coloc.tsv`
(**Supplementary Table S4**), `tox3_crossancestry_local_rg.csv` (Figure 5), `xanc_coloc_afr.csv`
(African-ancestry check at the colocalizing loci), `lava_pd_rls_afr.csv` and `lava_pd_rls_amr.csv`
(the African- and Admixed-American-ancestry RLS x PD scans), `known_vs_novel_blocks.csv`
(genome-wide-significance status of every screened block, the GWS column of Supplementary Table S3),
and `functional_annotation.tsv` (**Table 2**).

## Software

- **LDSC** v1.0.1 (Python-3 port), conda env: Python 3.9, numpy 1.23.5, scipy 1.9.3, pandas 1.4.4
- **LAVA** v0.1.5 - R 4.3.3
- **coloc** v5.2.3 (`coloc.abf` and `coloc.susie`) - R 4.3.3
- **susieR** v0.14.2 - R 4.3.3 (multi-signal fine-mapping behind coloc-SuSiE)
- **PLINK** v1.9 / v2.0
- **Open Targets Platform** release 26.09 (accessed November 2026) - functional annotation (Table 2)
- LD reference: 1000 Genomes Phase 3 (EUR / AFR / AMR); HapMap3 SNP list for munging

Exact package pins and the LDSC Python-3 patches are recorded in `scripts/1_setup/00_setup_environment.sh`.

## Configuration

Every script refers to the analysis root through the placeholder **`/path/to/analysis`** (and the
private RLS input through `/path/to/akcimen_rls_sumstats`). Point them at your own directories
before running, for example:

    grep -rl /path/to/analysis scripts/ run_all.sh | xargs sed -i "s#/path/to/analysis#$PWD/work#g"

Under the analysis root the scripts expect `data/gwas/`, `data/ldsc_ref/`, `data/reference/`,
`results/`, and `tools/ldsc/` (the stage-1 scripts create the reference and tool directories).

## Data

All GWAS summary statistics, with releases and accessions, are listed in **Supplementary Table S1
of the manuscript**. Every dataset is publicly available except the multi-ancestry RLS GWAS, which was
provided by Akcimen et al. Raw sumstats and LD panels are omitted here for size and licensing; the
download scripts record their sources.

## Reproducibility

These are the analysis scripts as run, in a Linux (WSL Ubuntu) environment. The public-data steps
can be reproduced after installing the toolchain (stage 1) and repointing the paths; the
RLS-anchored steps additionally require the multi-ancestry RLS summary statistics, which are
access-controlled through Akcimen et al., so full end-to-end re-execution depends on obtaining that
data. Every step, parameter, and data source is documented for transparency.

## License

MIT - see [LICENSE](LICENSE).
