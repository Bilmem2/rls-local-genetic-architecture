#!/usr/bin/env bash
# =============================================================================
# pd_02 — Munge PD GWAS for RLS×PD (Nalls discovery + FinnGen replication). Short; AI runs.
#   Nalls harmonised: rsid, beta = real log-OR (+SE), per-SNP N_cases/N_controls -> --N-cas-col/--N-con-col.
#   FinnGen R11: rsids (may be multi -> first), alt=effect, ref=other, beta(+sebeta), pval; fixed N.
# =============================================================================
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
ROOT=/path/to/analysis
LDSC="$ROOT/tools/ldsc"
HM3="$ROOT/data/ldsc_ref/w_hm3.snplist"
PD="$ROOT/data/gwas/pd"
M="$ROOT/results/munged"; mkdir -p "$M"

# ---- FinnGen R11 G6_PARKINSON case/control N — FILL from pd_01 manifest output / risteys.finngen.fi ----
FG_NCAS=5150       # FinnGen R11 G6_PARKINSON num_cases  (manifest)
FG_NCON=448583     # FinnGen R11 G6_PARKINSON num_controls (manifest)
# ------------------------------------------------------------------------------------------------------

echo "### pd_02_munge — $(date -u +%FT%TZ) ###"

echo "=== 1. PD Nalls 2019 (discovery) ==="
if [ -s "$PD/nalls2019_PD_eur.h.tsv.gz" ]; then
  python "$LDSC/munge_sumstats.py" --sumstats "$PD/nalls2019_PD_eur.h.tsv.gz" \
    --snp rsid --a1 effect_allele --a2 other_allele --signed-sumstats beta,0 --p p_value \
    --N-cas-col N_cases --N-con-col N_controls \
    --merge-alleles "$HM3" --chunksize 500000 --out "$M/pd_nalls" > "$M/pd_nalls.mungelog" 2>&1 \
    && grep -E "remain|Median|Writing|Mean chi|WARNING: median" "$M/pd_nalls.mungelog" | tail -8 \
    || { echo "  !! FAIL"; tail -25 "$M/pd_nalls.mungelog"; }
else echo "  !! missing $PD/nalls2019_PD_eur.h.tsv.gz — run pd_01 first"; fi

echo "=== 2. PD FinnGen R11 (replication) ==="
if [ ! -s "$PD/finngen_R11_PD.gz" ]; then
  echo "  !! missing FinnGen file — run pd_01 first"
elif [ "$FG_NCAS" -gt 0 ] && [ "$FG_NCON" -gt 0 ]; then
  CLEAN="$PD/finngen_R11_PD.clean.tsv.gz"
  if [ ! -s "$CLEAN" ]; then
    echo "  cleaning multi-rsid + header # ..."
    zcat "$PD/finngen_R11_PD.gz" \
      | awk 'BEGIN{FS=OFS="\t"} NR==1{sub(/^#/,"",$1); print; next} {n=split($5,a,","); $5=a[1]; print}' \
      | gzip > "$CLEAN"
  fi
  python "$LDSC/munge_sumstats.py" --sumstats "$CLEAN" \
    --snp rsids --a1 alt --a2 ref --signed-sumstats beta,0 --p pval \
    --N-cas "$FG_NCAS" --N-con "$FG_NCON" \
    --merge-alleles "$HM3" --chunksize 500000 --out "$M/pd_finngen" > "$M/pd_finngen.mungelog" 2>&1 \
    && grep -E "remain|Median|Writing|Mean chi|WARNING: median" "$M/pd_finngen.mungelog" | tail -8 \
    || { echo "  !! FAIL"; tail -25 "$M/pd_finngen.mungelog"; }
else
  echo "  !! FG_NCAS/FG_NCON still 0 — fill them (pd_01 prints the FinnGen N), then rerun pd_02."
fi

echo "=== munged PD outputs ==="
ls -la "$M"/pd_*.sumstats.gz 2>/dev/null
echo "=== PD_MUNGE_DONE ==="
