#!/usr/bin/env bash
# ============================================================================
# Canonical end-to-end execution order for the RLS local-genetic-architecture
# pipeline. Set the paths first (README > Configuration): every script uses the
# placeholder /path/to/analysis.
#
# These are the analysis scripts as run. The public-data steps are reproducible;
# the RLS-anchored steps additionally require the access-controlled multi-ancestry
# RLS summary statistics (Akcimen et al.). This driver documents the order; remove
# the leading "#" on a step to skip it, or run individual scripts directly.
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
step(){ echo; echo ">>> $*"; "$@"; }

# ---- Stage 1: environment + LD references ----
step bash    scripts/1_setup/00_setup_environment.sh
step bash    scripts/1_setup/01b_get_ldsc_ref.sh
step bash    scripts/1_setup/01c_smoketest.sh
step bash    scripts/1_setup/04_get_1000g_eur.sh
step bash    scripts/1_setup/dl_g1000_panels.sh
step bash    scripts/1_setup/build_hm3_hg38_map.sh

# ---- Stage 2: download + munge GWAS ----
step bash    scripts/2_munge/16_munge_akcimen.sh          # RLS EUR
step bash    scripts/2_munge/dl_rls_xanc.sh               # RLS multi-ancestry (download)
step bash    scripts/2_munge/munge_rls_xanc.sh            # RLS multi-ancestry (munge)
step bash    scripts/2_munge/pd_01_download.sh            # PD Nalls + FinnGen (download)
step bash    scripts/2_munge/pd_02_munge.sh               # PD Nalls (munge)
step bash    scripts/2_munge/munge_pd_rizig.sh            # PD AFR
step python  scripts/2_munge/annotate_pd_xanc.py          # annotate rsIDs on hg38 PD
step bash    scripts/2_munge/munge_pd_xanc.sh             # PD GP2 (EUR) + Loesch (AMR)
step bash    scripts/2_munge/munge_pd_largepd.sh          # PD LARGE-PD (AMR)
step bash    scripts/2_munge/munge_panel.sh               # psychiatric + sleep
step bash    scripts/2_munge/munge_insomnia.sh
step bash    scripts/2_munge/annotate_sleep.sh            # annotate rsIDs (OSA, sleep duration)
step bash    scripts/2_munge/munge_sleep2.sh              # OSA + sleep duration
step bash    scripts/2_munge/munge_sleep3.sh              # sleep-duration rsID fix
step bash    scripts/2_munge/munge_rg_daytime.sh          # daytime sleepiness
step bash    scripts/2_munge/munge_narcolepsy_ollila.sh   # narcolepsy (Ollila 2023 EUR)
step bash    scripts/2_munge/dl_finngen_osa.sh            # FinnGen R12 obstructive sleep apnea (download)

# ---- Stage 3: global genetic correlation (LDSC) ----
step bash    scripts/3_global_rg/panel_rg.sh
step bash    scripts/3_global_rg/panel_rg_akcimen.sh      # primary (RLS-Akcimen as phenotype 1)
step bash    scripts/3_global_rg/insomnia_rg.sh
step bash    scripts/3_global_rg/pd_03_rg.sh

# ---- Stage 4: local genetic correlation (LAVA) ----
step bash    scripts/4_local_rg_lava/panel_lava_run.sh    # genome-wide (panel_lava2_run.sh = parallel alternative)
step Rscript scripts/4_local_rg_lava/run_narco_lava.R     # narcolepsy (Ollila 2023) against RLS and PD
step Rscript scripts/4_local_rg_lava/lava_test_family.R   # unique test family, Bonferroni + BH -> results/lava_test_family.tsv
step bash    scripts/4_local_rg_lava/pd_xanc_run_lava.sh  # Admixed-American RLS x PD scan
step bash    scripts/4_local_rg_lava/neff_from_se.sh      # effective N of the GP2 meta-analysis from its standard errors
step bash    scripts/4_local_rg_lava/xanc_lava_tox3.sh    # African-ancestry RLS x PD scan + GP2 at TOX3

# ---- Stage 5: colocalization ----
step Rscript scripts/5_coloc/coloc_screen_full.R               # coloc.abf over every block with LAVA P < 1e-3
step Rscript scripts/4_local_rg_lava/known_vs_novel_blocks.R   # GWS status of each screened block
step Rscript scripts/4_local_rg_lava/gws_full_files.R          # same, short/long sleep and OSA from the full files
step python  scripts/5_coloc/build_table_s3.py                 # Supplementary Table S3
step Rscript scripts/5_coloc/coloc_sensitivity.R               # window x prior + coloc-SuSiE -> Supp. Table S4, panels A-C
step Rscript scripts/5_coloc/finngen_osa_coloc.R               # RLS x FinnGen OSA -> Supp. Table S4, panel E

# ---- Stage 6: TOX3/CASC16 cross-ancestry ----
step Rscript scripts/6_cross_ancestry_tox3/run_xanc_coloc_afr.R  # AFR local rg at the colocalizing loci
step bash    scripts/6_cross_ancestry_tox3/betabeta_extract.sh   # also renders betabeta_plot.py
step python  scripts/6_cross_ancestry_tox3/forest_tox3.py       # Figure 5

# ---- Stage 7: figures ----
step Rscript scripts/7_figures/forest_rg.R                      # Figure 2
step python  scripts/7_figures/coloc_selectivity.py             # Figure 3

echo; echo "=== pipeline complete ==="
