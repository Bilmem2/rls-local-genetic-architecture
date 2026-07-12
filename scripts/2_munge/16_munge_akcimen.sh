#!/usr/bin/env bash
# Akçimen 2026 EUR RLS meta (36,993 ca / 639,182 co, eff N~139,900) — munge by rsID.
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
ROOT=/path/to/analysis
LDSC=$ROOT/tools/ldsc
HM3=$ROOT/data/ldsc_ref/w_hm3.snplist
M=$ROOT/results/munged
SRC="/path/to/akcimen_rls_sumstats"
DST=$ROOT/data/gwas/rls/akcimen2026_EUR_RLS.txt

echo "=== copy to native fs ==="
[ -f "$DST" ] || cp "$SRC" "$DST"
ls -lh "$DST"; echo "header:"; head -1 "$DST"

echo "=== munge (rsID; cols SNPID CHR POS EA NEA BETA SE P STATUS rsID) ==="
python "$LDSC/munge_sumstats.py" --sumstats "$DST" --N 139900 \
  --snp rsID --a1 EA --a2 NEA --p P --signed-sumstats BETA,0 \
  --ignore SNPID,STATUS,POS,CHR \
  --merge-alleles "$HM3" --chunksize 500000 \
  --out "$M/rls_akcimen" > "$M/rls_akcimen.mungelog" 2>&1 \
  && grep -E "remain|Median value|Writing summary|WARNING: median" "$M/rls_akcimen.mungelog" | tail -6 \
  || { echo "FAIL"; tail -25 "$M/rls_akcimen.mungelog"; }
echo "=== munged ==="; ls -la "$M/rls_akcimen.sumstats.gz" 2>/dev/null
echo "=== AKCIMEN_MUNGE_DONE ==="
