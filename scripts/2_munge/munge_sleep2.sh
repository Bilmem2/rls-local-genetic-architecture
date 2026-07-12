#!/usr/bin/env bash
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
cd /path/to/analysis || exit 1
OUT=results/munged; MA=data/ldsc_ref/w_hm3.snplist
SL=data/gwas/sleep; GEL=$SL/sleepdur_gelernter
LOG=$OUT/munge_sleep2.log; : > "$LOG"
echo "START $(date)" >> "$LOG"
M(){ echo "===== $1 =====" >>"$LOG"; local out=$2; shift 2; \
  python tools/ldsc/munge_sumstats.py --merge-alleles "$MA" --out "$out" "$@" >>"$LOG" 2>&1 \
  && echo "OK $out" >>"$LOG" || echo "FAIL $out" >>"$LOG"; }

# OSA (trans-ancestry, hg38-annotated): beta null=0, total N=980757
M osa "$OUT/osa" --sumstats "$SL/osa/osa_rsid.tsv.gz" --snp SNP --a1 effect_allele --a2 other_allele --p p_value --signed-sumstats beta,0 --N 980757

# Gelernter sleep duration (hg19-annotated): Effect=beta null=0, per-SNP N
for tag in short long; do for anc in eur afr; do
  M "sleepdur_${anc}_${tag}" "$OUT/sleepdur_${anc}_${tag}" --sumstats "$GEL/${anc}_${tag}_rsid.txt.gz" --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats Effect,0 --N-col N
done; done

echo "DONE $(date)" >> "$LOG"
grep -E '^(=====|OK |FAIL )' "$LOG"
echo "--- metrics ---"; grep -E "Writing summary|Mean chi" "$LOG"
