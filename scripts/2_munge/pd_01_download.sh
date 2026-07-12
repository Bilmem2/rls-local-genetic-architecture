#!/usr/bin/env bash
# =============================================================================
# pd_01 — Download PD GWAS for the RLS×PD analysis (Fulya collaboration).
#   >>> RUN IN YOUR OWN WSL TERMINAL <<<  (~1.1 GB; idempotent + resumable).
#
#   Discovery   : Nalls et al. 2019 EUR (GCST009325, +UKB proxy, excl-23andMe), GWAS Catalog
#                 harmonised file (has rsID; beta=real log-OR + SE; per-SNP N_cases/N_controls).
#   Replication : FinnGen R11 G6_PARKINSON (independent EUR/Finnish; NOT in IPDGC/Nalls).
#   Mirrors the RLS 3-GWAS design: discovery + independent replication for any local (LAVA) hit.
# =============================================================================
set -uo pipefail
ROOT=/path/to/analysis
PD="$ROOT/data/gwas/pd"; mkdir -p "$PD"
LOG="$ROOT/results/_logs"; mkdir -p "$LOG"
exec > >(tee "$LOG/pd_01_download.log") 2>&1
echo "### pd_01_download — $(date -u +%FT%TZ) ###"

dl () {  # dl <url> <out> <expected_bytes|0>
  local url="$1" out="$2" want="${3:-0}"
  if [ -s "$out" ] && [ "$want" != "0" ]; then
    local h; h=$(stat -c%s "$out" 2>/dev/null || echo 0)
    [ "$h" = "$want" ] && { echo "  [skip,complete] $(basename "$out") ($h B)"; return 0; }
  fi
  echo "  [get] $(basename "$out")"
  curl -fL -C - --retry 5 --retry-delay 5 -o "$out" "$url" || { echo "  !! FAIL $url"; return 1; }
  echo "  [ok] $(stat -c%s "$out" 2>/dev/null) B"
}

echo "=== Discovery: Nalls 2019 EUR (GCST009325 harmonised) ==="
dl "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/harmonised/GCST009325.h.tsv.gz" \
   "$PD/nalls2019_PD_eur.h.tsv.gz" 0

echo "=== Replication: FinnGen R11 PD (G6_PARKINSON) ==="
dl "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_G6_PARKINSON.gz" \
   "$PD/finngen_R11_PD.gz" 812752514

echo "=== FinnGen R11 PD case/control N (needed for munge) ==="
for M in \
  "https://storage.googleapis.com/finngen-public-data-r11/finngen_R11_manifest.tsv" \
  "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_manifest.tsv" ; do
  if curl -sL "$M" 2>/dev/null | head -1 | grep -qi 'phenocode'; then
    echo "  manifest: $M"
    curl -sL "$M" 2>/dev/null | awk 'NR==1 || /G6_PARKINSON/' | grep -iE 'phenocode|G6_PARKINSON' | head -3
    break
  fi
done
echo "  (fallback: risteys.finngen.fi/endpoints/G6_PARKINSON shows n_cases / n_controls)"

echo "=== integrity + headers ==="
gzip -t "$PD/nalls2019_PD_eur.h.tsv.gz" 2>/dev/null && echo "  [ok] nalls gz" || echo "  !! nalls CORRUPT"
gzip -t "$PD/finngen_R11_PD.gz" 2>/dev/null && echo "  [ok] finngen gz" || echo "  !! finngen CORRUPT"
echo "  nalls header:";   zcat "$PD/nalls2019_PD_eur.h.tsv.gz" | head -1
echo "  finngen header:"; zcat "$PD/finngen_R11_PD.gz" | head -1
echo "### pd_01 DONE — next: pd_02 munge (AI runs) ###"
