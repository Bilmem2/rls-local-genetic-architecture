#!/usr/bin/env bash
# Munge Rizig 2023 African/admixed PD (has rsid -> build-agnostic, matches hg19 g1000_afr panel by rsID).
# N: 1,488 cases / 196,430 controls. Cols: chromosome base_pair_location effect_allele other_allele beta SE eaf p_value variant_id ref_allele direction ... rsid
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"; conda activate ldsc
ROOT=/path/to/analysis
LDSC=$ROOT/tools/ldsc; HM3=$ROOT/data/ldsc_ref/w_hm3.snplist
PD=$ROOT/data/gwas/pd; M=$ROOT/results/munged
RZ="$PD/rizig2023_PD_afr.tsv.gz"
echo "### extract Rizig from zip ###"
[ -s "$RZ" ] || unzip -p "$PD/Rizig_Parkinsons_2023.zip" "release5_11052023_summary_statistics_Rizig_et_al_2023_AFR_AAC_metaGWAS_no23andMe_hg38.txt.gz" > "$RZ"
echo "header: $(zcat "$RZ" | head -1)"
echo "### munge Rizig AFR (1488ca/196430co) ###"
python "$LDSC/munge_sumstats.py" --sumstats "$RZ" --N-cas 1488 --N-con 196430 \
  --snp rsid --a1 effect_allele --a2 other_allele --signed-sumstats beta,0 --p p_value \
  --merge-alleles "$HM3" --chunksize 500000 --out "$M/pd_rizig_afr" > "$M/pd_rizig_afr.mungelog" 2>&1 \
  && grep -E "remain|Median value|Writing summary|WARNING: median" "$M/pd_rizig_afr.mungelog" | tail -6 \
  || { echo "  !! FAIL"; tail -25 "$M/pd_rizig_afr.mungelog"; }
ls -la "$M/pd_rizig_afr.sumstats.gz" 2>/dev/null
echo "### RIZIG_MUNGE_DONE ###"
