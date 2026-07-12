#!/usr/bin/env bash
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
cd /path/to/analysis || exit 1
OUT=results/munged; MA=data/ldsc_ref/w_hm3.snplist
F=data/gwas/sleep/insomnia_watanabe/Watanabe_2022/insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz
LOG=$OUT/munge_insomnia.log; : > "$LOG"
python tools/ldsc/munge_sumstats.py --sumstats "$F" \
  --snp RSID_UKB --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 \
  --N-cas 109548 --N-con 277440 --ignore SNP,SNPID_UKB --merge-alleles "$MA" \
  --out "$OUT/insomnia" >> "$LOG" 2>&1
echo "EXIT $?" >> "$LOG"
grep -E "Writing summary|Mean chi|ERROR|expected to be|Found .* different|EXIT" "$LOG"
