#!/usr/bin/env bash
# Munge LARGE-PD 2025 (Leal) SAIGE P1+P2 Latino PD meta for AMR (replaces Loesch 2021).
# rs_number(rsID) + reference_allele(=effect) + other_allele + OR + p-value + per-SNP n_samples.
# Pre-filter to rsID rows with numeric OR/P/N (GWAMA output has some empty fields) -> then munge.
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"; conda activate ldsc
ROOT=/path/to/analysis; LDSC=$ROOT/tools/ldsc; HM3=$ROOT/data/ldsc_ref/w_hm3.snplist; M=$ROOT/results/munged
SRC="$ROOT/data/gwas/pd/LARGE-PD_PAPER/Meta/GWAMA/SAIGE_Meta-analysis_random_LARGE-PD_P1_P2_hg38.txt.gz"
CLEAN="$M/pd_largepd_amr.clean.tsv"

echo "### pre-filter: rsID + numeric OR(7)/p(12)/n_samples(18) ###"
printf 'SNP\tA1\tA2\tOR\tP\tN\n' > "$CLEAN"
zcat "$SRC" | awk -F'\t' 'NR>1 && $1 ~ /^rs/ && $7 ~ /^[0-9.eE+-]+$/ && $12 ~ /^[0-9.eE+-]+$/ && $18 ~ /^[0-9]+$/ && $18+0>0 {print $1"\t"$4"\t"$5"\t"$7"\t"$12"\t"$18}' >> "$CLEAN"
echo "kept rows: $(($(wc -l < "$CLEAN")-1))"
MEDN=$(awk 'NR>1{print $6}' "$CLEAN" | sort -n | awk '{a[NR]=$1} END{print (NR%2)?a[(NR+1)/2]:int((a[int(NR/2)]+a[int(NR/2)+1])/2)}')
echo "median N = $MEDN"; echo "$MEDN" > "$M/pd_largepd_amr.medN"

echo "### munge ###"
python "$LDSC/munge_sumstats.py" --sumstats "$CLEAN" --snp SNP --a1 A1 --a2 A2 --signed-sumstats OR,1 --p P --N-col N \
  --merge-alleles "$HM3" --chunksize 500000 --out "$M/pd_largepd_amr" > "$M/pd_largepd_amr.mungelog" 2>&1 \
  && grep -E "remain|Median value|Writing summary|WARNING: median|Mean chi" "$M/pd_largepd_amr.mungelog" | tail -8 \
  || { echo "  !! FAIL"; tail -25 "$M/pd_largepd_amr.mungelog"; }
ls -la "$M/pd_largepd_amr.sumstats.gz" 2>/dev/null
echo "### LARGEPD_MUNGE_DONE ###"
