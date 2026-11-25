#!/usr/bin/env bash
# FinnGen R12 obstructive sleep apnea (G6_SLEEPAPNO) summary statistics, used by
# scripts/5_coloc/finngen_osa_coloc.R. Source: FinnGen public data release R12.
set -eu
ROOT=/path/to/analysis
D=$ROOT/data/gwas/sleep/finngen_osa_r12
U=https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_G6_SLEEPAPNO.gz
mkdir -p "$D"; cd "$D"
curl -sS -L -m 3600 -o finngen_R12_G6_SLEEPAPNO.gz "$U"
gzip -t finngen_R12_G6_SLEEPAPNO.gz && echo GZIP_OK
zcat finngen_R12_G6_SLEEPAPNO.gz | awk 'NR<=2'
