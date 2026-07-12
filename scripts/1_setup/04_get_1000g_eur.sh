#!/usr/bin/env bash
# Phase 0/1/3: 1000G Phase3 EUR PLINK panel (MAGMA g1000_eur, hg19) for LD-clumping + LAVA.
set -uo pipefail
REFD=/path/to/analysis/data/reference
mkdir -p "$REFD"; cd "$REFD"

URL="https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download"
if [ ! -f g1000_eur.bed ]; then
  echo "downloading g1000_eur.zip (~488M)..."
  wget -q "$URL" -O g1000_eur.zip
  echo "unzipping..."
  unzip -o g1000_eur.zip >/dev/null 2>&1
  # files may extract flat or in a subdir; normalize
  if [ ! -f g1000_eur.bed ]; then
    found=$(find . -name 'g1000_eur.bed' | head -1)
    [ -n "$found" ] && mv "$(dirname "$found")"/g1000_eur.* . 2>/dev/null || true
  fi
fi

echo "=== contents ==="
ls -la g1000_eur.* 2>/dev/null
echo "=== variants (bim) ==="; wc -l < g1000_eur.bim 2>/dev/null
echo "=== samples (fam) ==="; wc -l < g1000_eur.fam 2>/dev/null
echo "=== bim head (chr snp cM bp a1 a2) ==="; head -3 g1000_eur.bim 2>/dev/null
echo "=== md5 of bed/bim/fam ==="; md5sum g1000_eur.bed g1000_eur.bim g1000_eur.fam 2>/dev/null
# tidy: keep zip removed to save space once extracted ok
[ -f g1000_eur.bed ] && rm -f g1000_eur.zip
echo "=== 1000G_DONE ==="
