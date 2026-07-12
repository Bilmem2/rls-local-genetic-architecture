#!/usr/bin/env bash
# rls_akcimen (primary RLS, biggest) x full panel — better-powered global rg + LAVA overlap intercepts
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
ROOT=/path/to/analysis
LDSC=$ROOT/tools/ldsc; REF=$ROOT/data/ldsc_ref; M=$ROOT/results/munged
OUT=$ROOT/results/ldsc_rg_panel; mkdir -p "$OUT"
REFLD="$REF/LDscore/LDscore."
WLD="$REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
PANEL="pd_nalls anx mdd bip ptsd scz rbd chronotype osa sleepdur_eur_short sleepdur_eur_long narcolepsy insomnia daytime_sleepiness"
blist(){ local l="$M/$1.sumstats.gz"; shift; for t in "$@"; do l="$l,$M/$t.sumstats.gz"; done; printf '%s' "$l"; }

echo "START $(date)"
python "$LDSC/ldsc.py" --rg "$(blist rls_akcimen $PANEL)" \
  --ref-ld-chr "$REFLD" --w-ld-chr "$WLD" --out "$OUT/rls_akcimen_vs_panel" > "$OUT/rls_akcimen_vs_panel.run.log" 2>&1 && echo "OK" || echo "FAIL"
echo "DONE $(date)"
echo "===== rls_akcimen own h2 / intercept (phenotype 1) ====="
grep -A7 "Heritability of phenotype 1$" "$OUT/rls_akcimen_vs_panel.log" 2>/dev/null | grep -E "Total Observed|Intercept|Lambda|Mean Chi"
echo "===== Summary of rg (h2_int = trait intercept; gcov_int = overlap) ====="
sed -n '/Summary of Genetic Correlation/,$p' "$OUT/rls_akcimen_vs_panel.log" 2>/dev/null
