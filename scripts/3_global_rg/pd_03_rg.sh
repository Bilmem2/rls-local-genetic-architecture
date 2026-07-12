#!/usr/bin/env bash
# =============================================================================
# pd_03 — Global LDSC genetic correlation: PD × RLS. Short; AI runs.
#   PD: Nalls 2019 (discovery) + FinnGen R11 (replication).
#   RLS: Akçimen meta (primary), Didriksen, Schormair.
#   Updates Estiar 2021 (rg=-0.028, null) with the larger RLS meta + PD discovery+replication.
#   Also prints the h2 intercepts + gencov intercept needed for the LAVA overlap matrix (pd_05).
# =============================================================================
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
ROOT=/path/to/analysis
LDSC="$ROOT/tools/ldsc"; REF="$ROOT/data/ldsc_ref"; M="$ROOT/results/munged"
OUT="$ROOT/results/ldsc_rg_pd"; mkdir -p "$OUT"
REFLD="$REF/LDscore/LDscore."
WLD="$REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."

rg () {  # <label> <pd_prefix> <rls_prefix>
  python "$LDSC/ldsc.py" --rg "$M/$2.sumstats.gz,$M/$3.sumstats.gz" \
    --ref-ld-chr "$REFLD" --w-ld-chr "$WLD" --out "$OUT/$1" > "$OUT/$1.run.log" 2>&1 \
    || { echo "  !! FAIL $1"; tail -15 "$OUT/$1.run.log"; }
}

echo "### pd_03_rg — $(date -u +%FT%TZ) ###"
LABELS=()
for pd in nalls finngen; do
  [ -s "$M/pd_${pd}.sumstats.gz" ] || { echo "  [skip] pd_${pd} not munged yet"; continue; }
  for rls in rls_akcimen rls_didriksen rls; do
    [ -s "$M/${rls}.sumstats.gz" ] || continue
    lab="${pd}_X_${rls}"
    echo "=== rg $lab ==="
    rg "$lab" "pd_${pd}" "$rls"
    LABELS+=("$lab")
  done
done

echo
echo "================= rg SUMMARY (PD × RLS) ================="
printf '%-26s %8s %7s %9s %11s\n' pair rg se p gcov_int
for l in "${LABELS[@]}"; do
  f="$OUT/$l.log"
  [ -f "$f" ] || { printf '%-26s   (no .log)\n' "$l"; continue; }
  read -r p1 p2 rg se z p h2 h2se h2int h2intse gcov gcovse <<< "$(tail -1 "$f")"
  printf '%-26s %8s %7s %9s %11s\n' "$l" "$rg" "$se" "$p" "$gcov"
done
echo
echo "--- h2 + intercepts per pair (for LAVA overlap matrix in pd_05) ---"
for l in "${LABELS[@]}"; do
  echo "== $l =="
  grep -E "Total Observed scale h2|Intercept:|Genetic Correlation:|Z-score:|^P:|Mean z1\*z2" "$OUT/$l.log" 2>/dev/null
done
echo "=== PD_RG_DONE  (next: pd_04 MR, pd_05 LAVA) ==="
