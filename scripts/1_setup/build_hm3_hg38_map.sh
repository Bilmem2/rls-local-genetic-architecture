#!/usr/bin/env bash
# Build a HM3 hg38 -> rsID map so we can annotate rsIDs onto hg38 PD sumstats that lack them
# (GP2, Loesch). HM3 is the SNP set munge keeps anyway -> smallest sufficient liftover (~1.2M).
# map cols: chr  pos_hg38  rsID  A1  A2   (A1/A2 from g1000_eur.bim; allele-set match downstream)
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"; conda activate ldsc
ROOT=/path/to/analysis
REF=$ROOT/data/reference; HM3=$ROOT/data/ldsc_ref/w_hm3.snplist
W=$REF/xanc_map; mkdir -p "$W"; cd "$W"

echo "### 1) CrossMap ###"
python -c "import CrossMap" 2>/dev/null || { echo "[install CrossMap]"; pip install --quiet CrossMap 2>&1 | tail -3; }
command -v CrossMap >/dev/null 2>&1 && CM="CrossMap" || CM="python -m CrossMap"

echo "### 2) hg19ToHg38 chain ###"
[ -s hg19ToHg38.over.chain.gz ] || wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
ls -lah hg19ToHg38.over.chain.gz

echo "### 3) HM3 SNPs from g1000_eur.bim as hg19 BED ###"
awk '$1 ~ /^rs/{print $1}' "$HM3" | sort -u > hm3.rs
# g1000_eur.bim: chr rsID cm pos A1 A2 ; keep HM3 ; BED0 = pos-1..pos ; name = rsID|A1|A2
awk 'NR==FNR{h[$1]=1; next} ($2 in h){print "chr"$1"\t"($4-1)"\t"$4"\t"$2"|"$5"|"$6}' hm3.rs "$REF/g1000_eur.bim" > hm3_hg19.bed
echo "HM3 in panel: $(wc -l < hm3_hg19.bed)"

echo "### 4) liftover hg19 -> hg38 ###"
$CM bed hg19ToHg38.over.chain.gz hm3_hg19.bed hm3_hg38.bed > crossmap.log 2>&1 || { echo "CrossMap FAIL"; tail -15 crossmap.log; exit 1; }

echo "### 5) build map (chr pos38 rsID A1 A2) ###"
awk 'BEGIN{OFS="\t"}{split($4,a,"|"); c=$1; sub(/^chr/,"",c); print c, $3, a[1], a[2], a[3]}' hm3_hg38.bed | sort -k1,1 -k2,2n > hm3_hg38_map.tsv
echo "map lines: $(wc -l < hm3_hg38_map.tsv)"; echo "head:"; head -3 hm3_hg38_map.tsv
echo "### HM3_HG38_MAP_DONE ###"
