#!/usr/bin/env bash
# LDSC global genetic correlation: RLS x panel, PD x panel (free intercept -> sample overlap handled)
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
ROOT=/path/to/analysis
LDSC=$ROOT/tools/ldsc; REF=$ROOT/data/ldsc_ref; M=$ROOT/results/munged
OUT=$ROOT/results/ldsc_rg_panel; mkdir -p "$OUT"
REFLD="$REF/LDscore/LDscore."
WLD="$REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
PANEL="anx mdd bip ptsd scz rbd chronotype osa sleepdur_eur_short sleepdur_eur_long narcolepsy"

blist(){ local l="$M/$1.sumstats.gz"; shift; for t in "$@"; do l="$l,$M/$t.sumstats.gz"; done; printf '%s' "$l"; }

echo "START $(date)"
# RLS vs {PD + panel}
python "$LDSC/ldsc.py" --rg "$(blist rls pd_nalls $PANEL)" \
  --ref-ld-chr "$REFLD" --w-ld-chr "$WLD" --out "$OUT/rls_vs_panel" > "$OUT/rls_vs_panel.run.log" 2>&1 && echo "OK rls" || echo "FAIL rls"
# PD vs {panel}
python "$LDSC/ldsc.py" --rg "$(blist pd_nalls $PANEL)" \
  --ref-ld-chr "$REFLD" --w-ld-chr "$WLD" --out "$OUT/pd_vs_panel" > "$OUT/pd_vs_panel.run.log" 2>&1 && echo "OK pd" || echo "FAIL pd"
echo "DONE $(date)"
echo "===== RLS vs panel ====="; sed -n '/Summary of Genetic Correlation/,$p' "$OUT/rls_vs_panel.log" 2>/dev/null
echo "===== PD  vs panel ====="; sed -n '/Summary of Genetic Correlation/,$p' "$OUT/pd_vs_panel.log" 2>/dev/null
