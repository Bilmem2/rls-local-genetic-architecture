#!/usr/bin/env bash
# Extract TOX3/CASC16 region SNPs from raw PD + RLS sumstats, then plot beta/beta.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate ldsc
GW=~/adhd-rls-sex-stratified/data/gwas
W=~/adhd-rls-sex-stratified/results/betabeta; mkdir -p "$W"; cd "$W"

echo "[extract] PD Nalls chr16:52.04-53.39Mb (GRCh37); RLS chr16 (rsID-merged downstream)"
# PD Nalls: chromosome=1 base_pair_location=2 effect_allele=3 other_allele=4 beta=5 ... p_value=8 rsid=9
zcat "$GW/pd/nalls2019_PD_eur.h.tsv.gz" \
  | awk 'NR==1{next} $1==16 && $2>=52040000 && $2<=53390000 {print $9,$3,$4,$5,$8}' > pd_nalls_tox3.txt
# RLS Schormair (GCST90399568): same cols as Nalls
awk 'NR==1{next} $1==16 {print $9,$3,$4,$5,$8}' "$GW/rls/GCST90399568.tsv" > rls_schormair_chr16.txt
# RLS Akcimen: SNPID=1 CHR=2 POS=3 EA=4 NEA=5 BETA=6 SE=7 P=8 STATUS=9 rsID=10
awk 'NR==1{next} $2==16 {print $10,$4,$5,$6,$8}' "$GW/rls/akcimen2026_EUR_RLS.txt" > rls_akcimen_chr16.txt
# RLS Didriksen: Chr=1 PosB38=2 ID=3 rsName=4 OA=5 EA=6 EAFrq=7 Cohorts=8 OR=9 P=10
zcat "$GW/rls/didriksen2020_RLS_meta.txt.gz" \
  | awk 'NR==1{next} {c=$1; sub(/^chr/,"",c)} c==16 {print $4,$6,$5,$9,$10}' > rls_didriksen_chr16.txt

wc -l pd_nalls_tox3.txt rls_*_chr16.txt
echo "[plot]"
python ~/adhd-rls-sex-stratified/scripts/betabeta_plot.py
