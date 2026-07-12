#!/usr/bin/env bash
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
cd /path/to/analysis || exit 1
N=data/gwas/sleep/narcolepsy_proxy
LOG=results/munged/narco_fix3.log; : > "$LOG"
RE='^[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$'
# minimal clean 5-col file: SNP(rsid col9) A1(effect col3) A2(other col4) OR(odds_ratio col5) P(p_value col8)
zcat "$N/GCST90477562.tsv.gz" | awk -F'\t' -v re="$RE" 'BEGIN{OFS="\t"}
  NR==1{print "SNP","A1","A2","OR","P"; next}
  ($5 ~ re && $8 ~ re && $9 ~ /^rs/){print $9,$3,$4,$5,$8}' | gzip > "$N/narco_min.tsv.gz"
echo "rows: $(zcat "$N/narco_min.tsv.gz"|wc -l)" >>"$LOG"
echo "head:" >>"$LOG"; zcat "$N/narco_min.tsv.gz"|head -3 >>"$LOG"
python tools/ldsc/munge_sumstats.py --sumstats "$N/narco_min.tsv.gz" \
  --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 \
  --N-cas 1079 --N-con 449612 --merge-alleles data/ldsc_ref/w_hm3.snplist \
  --out results/munged/narcolepsy >> "$LOG" 2>&1
echo "EXIT $?" >>"$LOG"
grep -E "rows:|Writing summary|Mean chi|ERROR|expected to be|EXIT" "$LOG"
