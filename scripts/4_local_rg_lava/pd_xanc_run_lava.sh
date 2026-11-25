#!/usr/bin/env bash
# =============================================================================
# Admixed-American RLS x PD LAVA scan: Loesch et al. 2021 Latin-American PD x the Admixed-American arm of
# the RLS meta-analysis (g1000_amr). Genome-wide (2,495 blocks; both-trait univariate gate p < 0.05/2495;
# Bonferroni + BH) with pd_05_lava.R. Identity sample-overlap matrix; EUR-hg19 blocks, rsID-matched.
# The African-ancestry scan and the GP2 estimates at TOX3 are in xanc_lava_tox3.sh.
# =============================================================================
set -uo pipefail
ROOT=/path/to/analysis
L="$ROOT/results/lava"; M="$ROOT/results/munged"; mkdir -p "$L" "$ROOT/results/_logs"
LOC="$ROOT/data/reference/lava_blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
RS=/usr/bin/Rscript
FILT='Extracting|Aligning|Negative variance|out of bounds|Removing|Reading|Processing input|treated as BINARY|^\|'
mk_info(){ printf 'phenotype cases controls filename\n%s %s %s %s\n%s %s %s %s\n' "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" > "$1"; }
mk_ovl(){ printf '%s %s\n%s 1 0\n%s 0 1\n' "$2" "$3" "$2" "$3" > "$1"; }
runscan(){ # 1 tag  2 ref-basename  3 LAVA_PAIRS  4 out.csv
  export LAVA_REF="$ROOT/data/reference/$2"; export LAVA_LOCF="$LOC"; export LAVA_OUT="$4"; export LAVA_PAIRS="$3"
  echo "### $1  (ref=$2)  -> $4 ###"; cd "$ROOT"
  stdbuf -oL -eL "$RS" scripts/4_local_rg_lava/pd_05_lava.R 2>&1 | tee "$ROOT/results/_logs/lava_$1.log" | grep --line-buffered -vE "$FILT"
}
echo "### ADMIXED-AMERICAN LAVA  $(date -u +%FT%TZ) ###"
mk_info "$L/info.amr.txt" pd_amr 740 741 "$M/pd_loesch_amr.sumstats.gz" rls_amr 2024 91902 "$M/rls_amr.sumstats.gz"; mk_ovl "$L/ovl.amr.txt" pd_amr rls_amr
runscan AMR g1000_amr "amr,$L/info.amr.txt,$L/ovl.amr.txt,pd_amr,rls_amr" "$L/lava_pd_rls_amr.csv"
echo "### done: results/lava/lava_pd_rls_amr.csv ###"
