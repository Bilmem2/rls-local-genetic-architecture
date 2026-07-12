#!/usr/bin/env bash
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
cd /path/to/analysis || exit 1
REF=data/reference
SL=data/gwas/sleep
GEL=$SL/sleepdur_gelernter
HM3=data/ldsc_ref/w_hm3.snplist
LOG=$SL/annotate.log; : > "$LOG"
echo "START $(date)" >> "$LOG"

# ---- HM3 rsID set (munge keeps only these) ----
awk '$1 ~ /^rs/{print $1}' "$HM3" | sort > /tmp/hm3.rs
echo "HM3 rsIDs: $(wc -l < /tmp/hm3.rs)" >> "$LOG"

# ---- HM3-restricted position->rsID maps: key=chr:pos  rsid A1 A2 ----
# hg19 EUR/AFR from g1000 bims (chr rsID cM pos A1 A2)
awk 'NR==FNR{h[$1]=1;next} ($2 in h){print $1":"$4"\t"$2"\t"$5"\t"$6}' /tmp/hm3.rs "$REF/g1000_eur.bim" > /tmp/hg19_eur.map
awk 'NR==FNR{h[$1]=1;next} ($2 in h){print $1":"$4"\t"$2"\t"$5"\t"$6}' /tmp/hm3.rs "$REF/g1000_afr.bim" > /tmp/hg19_afr.map
# hg38 (already HM3) from hm3_hg38_map.tsv (chr pos38 rsID A1 A2)
awk 'BEGIN{OFS="\t"}{print $1":"$2, $3, $4, $5}' "$REF/xanc_map/hm3_hg38_map.tsv" > /tmp/hg38.map
echo "maps: hg19_eur=$(wc -l </tmp/hg19_eur.map) hg19_afr=$(wc -l </tmp/hg19_afr.map) hg38=$(wc -l </tmp/hg38.map)" >> "$LOG"

reader(){ case "$1" in *.gz) zcat "$1";; *) cat "$1";; esac; }

# annotate: prepend rsID to allele-set-matched lines; report matched/total to stderr
# args: mapfile infile chrcol poscol a1col a2col  -> stdout
annotate(){
  awk -v cc="$3" -v pc="$4" -v ac="$5" -v oc="$6" 'BEGIN{FS=OFS="\t"}
    NR==FNR{ rs[$1]=$2; A[$1]=toupper($3); B[$1]=toupper($4); next }
    FNR==1{ print "SNP", $0; next }
    { tot++; k=$cc":"$pc; if(k in rs){ ea=toupper($ac); ob=toupper($oc);
        if((ea==A[k]&&ob==B[k])||(ea==B[k]&&ob==A[k])){ print rs[k], $0; hit++ } } }
    END{ printf "  matched %d / %d\n", hit+0, tot+0 > "/dev/stderr" }' \
    "$1" <(reader "$2")
}

# ===== OSA (expected hg38): cols chr(1) pos(2) EA(3) OA(4) beta se eaf p =====
echo "== OSA vs hg38 ==" >> "$LOG"
annotate /tmp/hg38.map "$SL/osa/GCST90693193.tsv" 1 2 3 4 2>>"$LOG" | gzip > "$SL/osa/osa_rsid.tsv.gz"

# ===== Gelernter build sanity on EUR-short: cols rsid CHR(2) BP(3) A1(4) A2(5) Freq1 N Effect P =====
echo "== Gelernter EUR-short vs hg19_eur (sanity) ==" >> "$LOG"
annotate /tmp/hg19_eur.map "$GEL/eur_meta_short_dbgap.txt.gz" 2 3 4 5 2>>"$LOG" >/dev/null
echo "== Gelernter EUR-short vs hg38 (sanity) ==" >> "$LOG"
annotate /tmp/hg38.map "$GEL/eur_meta_short_dbgap.txt.gz" 2 3 4 5 2>>"$LOG" >/dev/null

# ===== annotate all Gelernter with hg19 (EUR->eur bim, AFR->afr bim) =====
for tag in short long; do
  for anc in eur afr; do
    f="$GEL/${anc}_meta_${tag}_dbgap.txt.gz"
    [ -s "$f" ] || { echo "MISSING $f" >> "$LOG"; continue; }
    echo "== Gelernter ${anc} ${tag} vs hg19_${anc} ==" >> "$LOG"
    annotate /tmp/hg19_${anc}.map "$f" 2 3 4 5 2>>"$LOG" | gzip > "$GEL/${anc}_${tag}_rsid.txt.gz"
  done
done

echo "DONE $(date)" >> "$LOG"
echo "=== outputs ===" >> "$LOG"
ls -lh "$SL/osa/osa_rsid.tsv.gz" "$GEL"/*_rsid.txt.gz >> "$LOG" 2>&1
grep -E "matched|MISSING|maps:|^== " "$LOG"
