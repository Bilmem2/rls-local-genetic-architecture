#!/usr/bin/env bash
# Download Fulya's cross-ancestry RLS sumstats (AFR, AMR) from Zenodo record 20725243.
set -uo pipefail
RLS=/path/to/analysis/data/gwas/rls
cd "$RLS"
echo "### dl_rls_xanc $(date -u +%FT%TZ) ###"
echo "[dl] RLS_AFR_META (444 MB) ..."
wget -q -c "https://zenodo.org/records/20725243/files/RLS_AFR_META.tsv.gz?download=1" -O RLS_AFR_META.tsv.gz
echo "[dl] RLS_AMR_META (288 MB) ..."
wget -q -c "https://zenodo.org/records/20725243/files/RLS_AMR_META.tsv.gz?download=1" -O RLS_AMR_META.tsv.gz
echo "### sizes ###"; ls -lah RLS_AFR_META.tsv.gz RLS_AMR_META.tsv.gz
echo "### md5 (expect AFR d2943dd6..., AMR de71fbbd...) ###"; md5sum RLS_AFR_META.tsv.gz RLS_AMR_META.tsv.gz
echo "### DL_RLS_DONE ###"
