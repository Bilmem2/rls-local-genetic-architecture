#!/usr/bin/env bash
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
cd /path/to/analysis || exit 1
OUT=results/munged; MA=data/ldsc_ref/w_hm3.snplist
GEL=data/gwas/sleep/sleepdur_gelernter
LOG=$OUT/munge_sleep3.log; : > "$LOG"
echo "START $(date)" >> "$LOG"
M(){ echo "===== $1 =====" >>"$LOG"; local out=$2; shift 2; \
  python tools/ldsc/munge_sumstats.py --merge-alleles "$MA" --out "$out" "$@" >>"$LOG" 2>&1 \
  && echo "OK $out" >>"$LOG" || echo "FAIL $out" >>"$LOG"; }
for tag in short long; do for anc in eur afr; do
  f="$GEL/${anc}_${tag}_rsid.txt.gz"; f2="$GEL/${anc}_${tag}_rsid2.txt.gz"
  # rename the stale chr:pos 'rsid' header token so munge uses only the real-rsID SNP column
  zcat "$f" | sed '1s/rsid/OLDID/' | gzip > "$f2"
  M "sleepdur_${anc}_${tag}" "$OUT/sleepdur_${anc}_${tag}" --sumstats "$f2" --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats Effect,0 --N-col N
done; done
echo "DONE $(date)" >> "$LOG"
grep -E '^(=====|OK |FAIL )' "$LOG"
echo "--- metrics ---"; grep -E "Writing summary|Mean chi" "$LOG"
