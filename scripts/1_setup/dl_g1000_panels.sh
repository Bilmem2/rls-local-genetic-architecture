#!/usr/bin/env bash
# Download MAGMA/LAVA 1000G AFR + AMR reference panels (Build37), to match g1000_eur for cross-ancestry LAVA.
set -uo pipefail
REF=/path/to/analysis/data/reference
cd "$REF"
echo "### dl_g1000_panels $(date -u +%FT%TZ) ###"
echo "[dl] g1000_afr.zip (1.0 GB) ..."
wget -q -c "https://vu.data.surfsara.nl/index.php/s/ePXET6IWVTwTes4/download" -O g1000_afr.zip
echo "[dl] g1000_amr.zip (492 MB) ..."
wget -q -c "https://vu.data.surfsara.nl/index.php/s/TXDEm70eEO7AgOb/download" -O g1000_amr.zip
echo "[unzip]"
unzip -o g1000_afr.zip >/dev/null 2>&1 && echo "  afr unzipped"
unzip -o g1000_amr.zip >/dev/null 2>&1 && echo "  amr unzipped"
echo "### resulting panel files ###"
ls -lah g1000_afr.* g1000_amr.* 2>/dev/null
echo "### DL_PANELS_DONE ###"
