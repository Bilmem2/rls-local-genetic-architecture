#!/usr/bin/env bash
HERE="$(cd "$(dirname "$0")" && pwd)"  # this script's directory
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
cd /path/to/analysis || exit 1
W=/path/to/analysis/coloc_rlspd; mkdir -p "$W"
RLS=data/gwas/rls/akcimen2026_EUR_RLS.txt
PD=data/gwas/pd/nalls2019_PD_eur.h.tsv.gz
# RLS akcimen: SNPID CHR POS EA NEA BETA SE P STATUS rsID  (tab)
# PD Nalls: chr pos ea oa beta se eaf p rsid Ncas Ncon  (tab)
extract(){ local c=$1 lo=$2 hi=$3 tag=$4
  awk -F'\t' -v c=$c -v lo=$lo -v hi=$hi '$2==c && $3>=lo && $3<=hi && $10~/^rs/ {print $10"\t"$4"\t"$5"\t"$6"\t"$7}' "$RLS" > "$W/rls_$tag.tsv"
  zcat "$PD" | awk -F'\t' -v c=$c -v lo=$lo -v hi=$hi 'NR>1 && $1==c && $2>=lo && $2<=hi && $9~/^rs/ {print $9"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > "$W/pd_$tag.tsv"
}
extract 11 133521643 134351064 chr11
extract 12 23923799  25058714  chr12
extract 16 52041337  53393882  tox3
echo "=== region SNP counts ==="; wc -l "$W"/*.tsv
Rscript "$HERE/coloc_rlspd.R"
