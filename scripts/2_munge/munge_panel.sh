#!/usr/bin/env bash
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
cd /path/to/analysis || exit 1
GW=data/gwas; PSY=$GW/psych; SL=$GW/sleep; OUT=results/munged
mkdir -p "$OUT"
MA=data/ldsc_ref/w_hm3.snplist
LOG=$OUT/munge_panel.log; : > "$LOG"
echo "START $(date)" >> "$LOG"

# preprocess VCF-style files: drop ## preamble, strip leading # from header
zcat "$PSY/eur_ptsd_pcs_v4_aug3_2021.vcf.gz" | grep -vE '^##' | sed '1s/^#//' | gzip > "$PSY/ptsd_eur_clean.tsv.gz"
zcat "$PSY/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz" | grep -vE '^##' | sed '1s/^#//' | gzip > "$PSY/scz_eur_clean.tsv.gz"

M() { echo "===== $1 =====" >> "$LOG"; local out="$2"; shift 2; \
  python tools/ldsc/munge_sumstats.py --merge-alleles "$MA" --out "$out" "$@" >> "$LOG" 2>&1 \
  && echo "OK $out" >> "$LOG" || echo "FAIL $out" >> "$LOG"; }

# psychiatric (DANER: OR null=1; per-SNP Nca/Nco)
M anxiety "$OUT/anx" --sumstats "$PSY/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz" --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --N-cas-col Nca --N-con-col Nco
M mdd "$OUT/mdd" --sumstats "$GW/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz" --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --N-cas-col Nca --N-con-col Nco
M bipolar "$OUT/bip" --sumstats "$PSY/bip2024_eur_noUKB_no23andMe.gz" --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats OR,1 --N-cas-col Nca --N-con-col Nco
# ptsd (Z null=0; NEFF)
M ptsd "$OUT/ptsd" --sumstats "$PSY/ptsd_eur_clean.tsv.gz" --snp ID --a1 A1 --a2 A2 --p P --signed-sumstats Z,0 --N-col NEFF
# scz (BETA null=0; NCAS/NCON)
M scz "$OUT/scz" --sumstats "$PSY/scz_eur_clean.tsv.gz" --snp ID --a1 A1 --a2 A2 --p PVAL --signed-sumstats BETA,0 --N-cas-col NCAS --N-con-col NCON
# sleep (rsID present)
M rbd "$OUT/rbd" --sumstats "$SL/rbd_krohn/GCST90204200_buildGRCh37.tsv" --snp variant_id --a1 effect_allele --a2 other_allele --p p_value --signed-sumstats beta,0 --N-col N
M chronotype "$OUT/chronotype" --sumstats "$SL/chronotype_jones/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz" --snp SNP --a1 ALLELE1 --a2 ALLELE0 --p P_BOLT_LMM --signed-sumstats BETA,0 --N 449734
M narcolepsy "$OUT/narcolepsy" --sumstats "$SL/narcolepsy_proxy/GCST90477562.tsv.gz" --snp rsid --a1 effect_allele --a2 other_allele --p p_value --signed-sumstats odds_ratio,1 --N-cas-col num_cases --N-con-col num_controls

echo "DONE $(date)" >> "$LOG"
echo "=== munged outputs ===" >> "$LOG"
ls -lh "$OUT"/*.sumstats.gz >> "$LOG" 2>&1
grep -E '^(=====|OK |FAIL )' "$LOG"
