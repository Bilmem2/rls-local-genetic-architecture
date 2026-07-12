#!/usr/bin/env bash
HERE="$(cd "$(dirname "$0")" && pwd)"  # this script's directory
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
cd /path/to/analysis || exit 1
W=/path/to/analysis/coloc_meis1; mkdir -p "$W"
RLS=data/gwas/rls/akcimen2026_EUR_RLS.txt
INS=data/gwas/sleep/insomnia_watanabe/Watanabe_2022/insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz
C=2; LO=65938003; HI=67224027
# RLS akcimen: SNPID CHR POS EA NEA BETA SE P STATUS rsID
awk -F'\t' -v c=$C -v lo=$LO -v hi=$HI '$2==c && $3>=lo && $3<=hi && $10~/^rs/ {print $10"\t"$4"\t"$5"\t"$6"\t"$7}' "$RLS" > "$W/rls.tsv"
# insomnia: SNP CHR BP A1 TEST NMISS OR SE ... A2(13) SNPID_UKB(14) RSID_UKB(15) ... MAF(20)
zcat "$INS" | awk -F'\t' -v c=$C -v lo=$LO -v hi=$HI 'NR>1 && $2==c && $3>=lo && $3<=hi && $15~/^rs/ && $7>0 {print $15"\t"$4"\t"$13"\t"log($7)"\t"$8"\t"$20}' > "$W/ins.tsv"
echo "=== SNP counts ==="; wc -l "$W"/*.tsv
Rscript "$HERE/coloc_meis1.R"
