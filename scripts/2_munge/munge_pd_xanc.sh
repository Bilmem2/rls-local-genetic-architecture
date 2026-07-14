#!/usr/bin/env bash
# Annotate (rsID via HM3 hg38 map) + munge GP2 EUR + Loesch AMR PD sumstats.
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"; conda activate ldsc
ROOT=/path/to/analysis; LDSC=$ROOT/tools/ldsc; HM3=$ROOT/data/ldsc_ref/w_hm3.snplist; M=$ROOT/results/munged

echo "### annotate rsIDs (hg38 map) ###"
python "$ROOT/scripts/2_munge/annotate_pd_xanc.py"
GP2N=$(awk '$1=="GP2"{print $2}' "$M/pd_xanc.N"); LON=$(awk '$1=="LOESCH"{print $2}' "$M/pd_xanc.N")
echo "derived N: GP2=$GP2N  Loesch=$LON"

echo "### munge GP2 EUR (N=$GP2N) ###"
python "$LDSC/munge_sumstats.py" --sumstats "$M/pd_gp2_eur.clean.tsv" --N "$GP2N" \
  --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P --merge-alleles "$HM3" --chunksize 500000 \
  --out "$M/pd_gp2_eur" > "$M/pd_gp2_eur.mungelog" 2>&1 \
  && grep -E "remain|Median value|Writing summary|WARNING: median" "$M/pd_gp2_eur.mungelog" | tail -6 \
  || { echo "  !! GP2 FAIL"; tail -25 "$M/pd_gp2_eur.mungelog"; }

echo "### munge Loesch AMR (N=$LON) ###"
python "$LDSC/munge_sumstats.py" --sumstats "$M/pd_loesch_amr.clean.tsv" --N "$LON" \
  --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P --merge-alleles "$HM3" --chunksize 500000 \
  --out "$M/pd_loesch_amr" > "$M/pd_loesch_amr.mungelog" 2>&1 \
  && grep -E "remain|Median value|Writing summary|WARNING: median" "$M/pd_loesch_amr.mungelog" | tail -6 \
  || { echo "  !! Loesch FAIL"; tail -25 "$M/pd_loesch_amr.mungelog"; }

echo "### munged ###"; ls -la "$M"/pd_gp2_eur.sumstats.gz "$M"/pd_loesch_amr.sumstats.gz 2>/dev/null
echo "### PD_XANC_MUNGE_DONE ###"
