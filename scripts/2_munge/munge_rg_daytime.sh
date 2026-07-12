#!/usr/bin/env bash
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
ROOT=/path/to/analysis
LDSC=$ROOT/tools/ldsc; REF=$ROOT/data/ldsc_ref; M=$ROOT/results/munged
RGOUT=$ROOT/results/ldsc_rg_panel; mkdir -p "$RGOUT"
MA=$REF/w_hm3.snplist
REFLD="$REF/LDscore/LDscore."
WLD="$REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
F=$ROOT/data/gwas/sleep/daytime_saxena/Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt
LOG=$M/munge_daytime.log; : > "$LOG"

# munge (BOLT-LMM: BETA null=0; UKB EUR N=452071)
python "$LDSC/munge_sumstats.py" --sumstats "$F" \
  --snp SNP --a1 ALLELE1 --a2 ALLELE0 --p P --signed-sumstats BETA,0 --N 452071 \
  --merge-alleles "$MA" --out "$M/daytime_sleepiness" >> "$LOG" 2>&1
echo "MUNGE EXIT $?" >> "$LOG"
grep -E "Writing summary|Mean chi|ERROR|expected to be|Found .* different|MUNGE EXIT" "$LOG"

# rg vs RLS and PD
if [ -s "$M/daytime_sleepiness.sumstats.gz" ]; then
  for anc in rls pd_nalls; do
    python "$LDSC/ldsc.py" --rg "$M/$anc.sumstats.gz,$M/daytime_sleepiness.sumstats.gz" \
      --ref-ld-chr "$REFLD" --w-ld-chr "$WLD" --out "$RGOUT/${anc}_vs_daytime" > "$RGOUT/${anc}_vs_daytime.run.log" 2>&1 && echo "rg OK $anc" || echo "rg FAIL $anc"
  done
  echo "=== RLS x daytime ==="; sed -n '/Summary of Genetic Correlation/,$p' "$RGOUT/rls_vs_daytime.log" 2>/dev/null | grep -Ei "daytime|p1 "
  echo "=== PD  x daytime ==="; sed -n '/Summary of Genetic Correlation/,$p' "$RGOUT/pd_nalls_vs_daytime.log" 2>/dev/null | grep -Ei "daytime|p1 "
fi
