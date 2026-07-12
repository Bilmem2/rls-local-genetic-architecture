#!/usr/bin/env bash
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
ROOT=/path/to/analysis
LDSC=$ROOT/tools/ldsc; REF=$ROOT/data/ldsc_ref; M=$ROOT/results/munged
OUT=$ROOT/results/ldsc_rg_panel; mkdir -p "$OUT"
REFLD="$REF/LDscore/LDscore."
WLD="$REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
for anc in rls pd_nalls; do
  python "$LDSC/ldsc.py" --rg "$M/$anc.sumstats.gz,$M/insomnia.sumstats.gz" \
    --ref-ld-chr "$REFLD" --w-ld-chr "$WLD" --out "$OUT/${anc}_vs_insomnia" > "$OUT/${anc}_vs_insomnia.run.log" 2>&1 \
    && echo "OK $anc" || echo "FAIL $anc"
done
echo "=== RLS x insomnia ==="; sed -n '/Summary of Genetic Correlation/,$p' "$OUT/rls_vs_insomnia.log" 2>/dev/null | grep -E "insomnia|rg|p1"
echo "=== PD x insomnia ==="; sed -n '/Summary of Genetic Correlation/,$p' "$OUT/pd_nalls_vs_insomnia.log" 2>/dev/null | grep -E "insomnia|rg|p1"
