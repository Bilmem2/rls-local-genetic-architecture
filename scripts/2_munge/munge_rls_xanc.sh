#!/usr/bin/env bash
# Munge cross-ancestry RLS (Akçimen multi-ancestry, medRxiv 2026.04.28). Same cols as EUR akcimen.
# Eff-N: AFR 2176ca/153313co -> 8581 ; AMR 2024ca/91902co -> 7922  (matches EUR --N 139900 convention)
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
ROOT=/path/to/analysis
LDSC=$ROOT/tools/ldsc
HM3=$ROOT/data/ldsc_ref/w_hm3.snplist
RLS=$ROOT/data/gwas/rls
M=$ROOT/results/munged

echo "### munge RLS AFR (eff N 8581) $(date -u +%FT%TZ) ###"
python "$LDSC/munge_sumstats.py" --sumstats "$RLS/RLS_AFR_META.tsv.gz" --N 8581 \
  --snp rsID --a1 EA --a2 NEA --p P --signed-sumstats BETA,0 --ignore SNPID,STATUS,POS,CHR \
  --merge-alleles "$HM3" --chunksize 500000 --out "$M/rls_afr" > "$M/rls_afr.mungelog" 2>&1 \
  && grep -E "remain|Median value|Writing summary|WARNING: median" "$M/rls_afr.mungelog" | tail -6 \
  || { echo "  !! AFR FAIL"; tail -25 "$M/rls_afr.mungelog"; }

echo "### munge RLS AMR (eff N 7922) ###"
python "$LDSC/munge_sumstats.py" --sumstats "$RLS/RLS_AMR_META.tsv.gz" --N 7922 \
  --snp rsID --a1 EA --a2 NEA --p P --signed-sumstats BETA,0 --ignore SNPID,STATUS,POS,CHR \
  --merge-alleles "$HM3" --chunksize 500000 --out "$M/rls_amr" > "$M/rls_amr.mungelog" 2>&1 \
  && grep -E "remain|Median value|Writing summary|WARNING: median" "$M/rls_amr.mungelog" | tail -6 \
  || { echo "  !! AMR FAIL"; tail -25 "$M/rls_amr.mungelog"; }

echo "### munged outputs ###"; ls -la "$M"/rls_afr.sumstats.gz "$M"/rls_amr.sumstats.gz 2>/dev/null
echo "### RLS_XANC_MUNGE_DONE ###"
