#!/usr/bin/env bash
ROOT=/path/to/analysis
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
LDSC=$ROOT/tools/ldsc; HM3=$ROOT/data/ldsc_ref/w_hm3.snplist
REF=$ROOT/data/ldsc_ref; M=$ROOT/results/munged; OUT=$ROOT/results/ldsc_rg_panel
REFLD="$REF/LDscore/LDscore."; WLD="$REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
ND=$ROOT/data/gwas/sleep/narcolepsy_ollila; mkdir -p "$ND"
zip=$ROOT/data/gwas/doi_10_5061_dryad_kd51c5b9b__v20221222.zip
F="CHR_Whites_info0.7.noHLA.noSig.hetP1e-5.LDSC.N.txt.gz"
[ -f "$ND/$F" ] || unzip -o "$zip" "$F" -d "$ND" >/dev/null

echo "=== munge Ollila EUR narcolepsy (beta,0; N-col Neff) ==="
python "$LDSC/munge_sumstats.py" --sumstats "$ND/$F" \
  --snp rsid --a1 allele_A --a2 allele_B --p P_value --signed-sumstats beta,0 \
  --N-col Neff --merge-alleles "$HM3" --chunksize 500000 \
  --out "$M/narcolepsy_ollila" > "$M/narcolepsy_ollila.mungelog" 2>&1 \
  && grep -E "SNPs remain|Median value|Writing summary|Mean chi|WARNING: median" "$M/narcolepsy_ollila.mungelog" | tail -6 \
  || { echo "MUNGE FAIL"; tail -20 "$M/narcolepsy_ollila.mungelog"; exit 1; }

for pair in "rls_akcimen" "pd_nalls"; do
  tag=$([ "$pair" = "pd_nalls" ] && echo "pd" || echo "rls_akcimen")
  python "$LDSC/ldsc.py" --rg "$M/$pair.sumstats.gz,$M/narcolepsy_ollila.sumstats.gz" \
    --ref-ld-chr "$REFLD" --w-ld-chr "$WLD" --out "$OUT/${tag}_vs_narco_ollila" >/dev/null 2>&1
done

echo; echo "=== NEW global rg with Ollila EUR narcolepsy ==="
for tag in rls_akcimen pd; do
  l="$OUT/${tag}_vs_narco_ollila.log"
  echo "-- ${tag} x narcolepsy_ollila --"
  awk '/Summary of Genetic Correlation/{f=1;next} f&&$1~/narcolepsy_ollila/{print "   rg="$3"  se="$4"  z="$5"  p="$6}' "$l"
  grep -m1 "Total Observed scale h2" "$l" | sed 's/^/   narco /'
done
echo "=== compare: OLD proxy narcolepsy rg was RLS 0.04 (p0.73), PD 0.01 (p0.93) ==="