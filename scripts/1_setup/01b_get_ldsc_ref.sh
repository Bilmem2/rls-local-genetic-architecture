#!/usr/bin/env bash
# Download current canonical LDSC reference (1000G Phase3 LD scores + HM3 weights + w_hm3 snplist)
set -uo pipefail
ROOT=/path/to/analysis
REF=$ROOT/data/ldsc_ref
mkdir -p "$REF"; cd "$REF"

echo "############ w_hm3.snplist ############"
if [ ! -f w_hm3.snplist ]; then
  wget -q "https://zenodo.org/records/7773502/files/w_hm3.snplist.gz?download=1" -O w_hm3.snplist.gz \
    && gunzip -f w_hm3.snplist.gz
fi
[ -f w_hm3.snplist ] && { echo "lines: $(wc -l < w_hm3.snplist)"; head -2 w_hm3.snplist; } || echo "FAILED w_hm3"

echo
echo "############ 1000G Phase3 LD scores (ref-ld) ############"
if [ ! -f 1000G_Phase3_ldscores.tgz ]; then
  wget -q "https://zenodo.org/records/10515792/files/1000G_Phase3_ldscores.tgz?download=1" -O 1000G_Phase3_ldscores.tgz
fi
tar -tzf 1000G_Phase3_ldscores.tgz >/dev/null 2>&1 && tar -xzf 1000G_Phase3_ldscores.tgz || echo "FAILED ldscores extract"

echo
echo "############ 1000G Phase3 HM3 weights (w-ld) ############"
if [ ! -f 1000G_Phase3_weights_hm3_no_MHC.tgz ]; then
  wget -q "https://zenodo.org/records/10515792/files/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1" -O 1000G_Phase3_weights_hm3_no_MHC.tgz
fi
tar -tzf 1000G_Phase3_weights_hm3_no_MHC.tgz >/dev/null 2>&1 && tar -xzf 1000G_Phase3_weights_hm3_no_MHC.tgz || echo "FAILED weights extract"

echo
echo "############ extracted structure ############"
find "$REF" -maxdepth 1 -mindepth 1 -type d -printf '%p\n'
echo "--- dirs with sample contents ---"
for d in "$REF"/*/; do
  n=$(find "$d" -type f | wc -l)
  echo "[$d]  files=$n"
  find "$d" -type f -printf '%f\n' | sort | head -4
done

echo
echo "############ md5 of downloads (for PROVENANCE) ############"
md5sum 1000G_Phase3_ldscores.tgz 1000G_Phase3_weights_hm3_no_MHC.tgz 2>/dev/null
echo "=== ALLDONE ==="
