#!/usr/bin/env bash
HERE="$(cd "$(dirname "$0")" && pwd)"  # this script's directory
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
cd /path/to/analysis || exit 1
W=/path/to/analysis/coloc_chr1; mkdir -p "$W"
G=data/gwas; C=1; LO=36855140; HI=38474036
# PD Nalls (tab): chr pos ea oa beta se eaf p rsid Ncas Ncon
zcat "$G/pd/nalls2019_PD_eur.h.tsv.gz" | awk -F'\t' -v c=$C -v lo=$LO -v hi=$HI 'BEGIN{OFS="\t";print "rsid","ea","oa","beta","se","eaf","Ncas","Ncon"} NR>1 && $1==c && $2>=lo && $2<=hi && $9~/^rs/ {print $9,$3,$4,$5,$6,$7,$10,$11}' > "$W/pd.tsv"
# MDD daner: CHR SNP BP A1 A2 FRQ_A FRQ_U INFO OR SE P
zcat "$G/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz" | awk -v c=$C -v lo=$LO -v hi=$HI 'BEGIN{OFS="\t";print "rsid","ea","oa","beta","se","eaf"} NR>1 && $1==c && $3>=lo && $3<=hi && $2~/^rs/ && $9>0 {print $2,$4,$5,log($9),$10,$7}' > "$W/mdd.tsv"
# anx daner
zcat "$G/psych/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz" | awk -v c=$C -v lo=$LO -v hi=$HI 'BEGIN{OFS="\t";print "rsid","ea","oa","beta","se","eaf"} NR>1 && $1==c && $3>=lo && $3<=hi && $2~/^rs/ && $9>0 {print $2,$4,$5,log($9),$10,$7}' > "$W/anx.tsv"
# ptsd vcf: #CHROM ID POS A1 A2 FREQ NEFF Z P
zcat "$G/psych/eur_ptsd_pcs_v4_aug3_2021.vcf.gz" | grep -vE '^##' | awk -F'\t' -v c=$C -v lo=$LO -v hi=$HI 'BEGIN{OFS="\t";print "rsid","ea","oa","z","freq","neff"} $1==c && $3>=lo && $3<=hi && $2~/^rs/ {print $2,$4,$5,$8,$6,$7}' > "$W/ptsd.tsv"
echo "=== region SNP counts (incl header) ==="; wc -l "$W"/*.tsv
Rscript "$HERE/coloc_chr1.R"
