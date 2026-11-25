#!/usr/bin/env bash
# =============================================================================
# African-ancestry RLS x PD LAVA scan and the GP2 estimates at TOX3/CASC16 (LAVA block 2134).
#   AFR: Rizig et al. 2023 PD, release without 23andMe: 1,200 cases / 2,445 controls (Table 1 of the
#        paper), x the African-ancestry arm of the RLS meta-analysis (g1000_afr). Genome-wide gated scan
#        with pd_05_lava.R; identity sample-overlap matrix (no shared cohorts).
#   GP2: European meta-analysis (Leonard et al. 2025; file GP2_ALL_EUR_ALL_DATASET_HG38_12162024) x the
#        three European RLS GWAS (g1000_eur). The file has no per-variant N; the effective N estimated
#        from the standard errors (neff_from_se.sh) is 170,212, entered as 85,106 cases / 85,106 controls.
#        GP2 and the RLS GWAS share biobank participants, so the overlap term is the LDSC cross-trait
#        intercept standardised as in the panel analysis: gcov_int / sqrt(h2_int_PD * h2_int_RLS).
# The N column of a munged file is constant across variants, so replacing it is equivalent to munging
# with that N (Z statistics are unchanged).
# =============================================================================
set -euo pipefail
ROOT=/path/to/analysis
REPO=$(cd "$(dirname "$0")/../.." && pwd)
M=$ROOT/results/munged; L=$ROOT/results/lava; mkdir -p "$L"
LOC=$ROOT/data/reference/lava_blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile
RS=/usr/bin/Rscript

echo "### 1. munged copies with the sample sizes above ###"
for f in pd_rizig_afr pd_gp2_eur; do
  h=$(zcat "$M/$f.sumstats.gz" | awk 'NR==1{print $5}'); [ "$h" = "N" ] || { echo "unexpected header in $f"; exit 1; }; done
zcat "$M/pd_rizig_afr.sumstats.gz" | awk 'BEGIN{FS=OFS="\t"} NR==1{print; next}{$5=3645; print}' | gzip > "$M/pd_rizig_afr_n3645.sumstats.gz"
zcat "$M/pd_gp2_eur.sumstats.gz"   | awk 'BEGIN{FS=OFS="\t"} NR==1{print; next}{$5=170212; print}' | gzip > "$M/pd_gp2_eur_n170212.sumstats.gz"

echo "### 2. LDSC: GP2 x the three European RLS GWAS (intercepts for the overlap term) ###"
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"; conda activate ldsc
python "$ROOT/tools/ldsc/ldsc.py" --rg "$M/pd_gp2_eur_n170212.sumstats.gz,$M/rls_akcimen.sumstats.gz,$M/rls_didriksen.sumstats.gz,$M/rls.sumstats.gz" \
  --ref-ld-chr "$ROOT/data/ldsc_ref/LDscore/LDscore." --w-ld-chr "$ROOT/data/ldsc_ref/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." \
  --out "$L/gp2_vs_rls" > /dev/null 2>&1
conda deactivate
# intercepts in log order: h2 GP2, then for each RLS GWAS: h2 RLS, gcov(GP2, RLS)
read -r H1 H2 G2 H3 G3 H4 G4 < <(awk '/^Intercept:/{printf "%s ", $2}' "$L/gp2_vs_rls.log" | awk '{print $1, $2, $3, $4, $5, $6, $7}')
std(){ awk -v g="$1" -v a="$2" -v b="$3" 'BEGIN{printf "%.6f", g/sqrt(a*b)}'; }
S_AKC=$(std "$G2" "$H1" "$H2"); S_DID=$(std "$G3" "$H1" "$H3"); S_SCH=$(std "$G4" "$H1" "$H4")
echo "standardised overlap: akcimen=$S_AKC didriksen=$S_DID schormair=$S_SCH"

mk_info(){ printf 'phenotype cases controls filename\n%s %s %s %s\n%s %s %s %s\n' "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" > "$1"; }
mk_ovl(){ printf '%s %s\n%s 1 %s\n%s %s 1\n' "$2" "$3" "$2" "$4" "$3" "$4" > "$1"; }
run(){ # 1 tag  2 reference  3 LAVA_PAIRS  4 output  5 locfile
  echo "### LAVA $1 (ref=$2, loci=$(($(wc -l < "$5")-1))) ###"
  LAVA_REF="$ROOT/data/reference/$2" LAVA_LOCF="$5" LAVA_OUT="$4" LAVA_PAIRS="$3" $RS "$REPO/scripts/4_local_rg_lava/pd_05_lava.R" > "$L/lava_$1.log" 2>&1 || true
  grep -E 'rho=|gate passed|input error' "$L/lava_$1.log" | awk 'NR<=20' || true; }

echo "### 3. GP2 at the TOX3 block ###"
head -1 "$LOC" > "$L/tox3.locfile"; awk '$1==2134' "$LOC" >> "$L/tox3.locfile"
mk_info "$L/info.gp2a.txt" pd_gp2 85106 85106 "$M/pd_gp2_eur_n170212.sumstats.gz" rls_akcimen 36993 639182 "$M/rls_akcimen.sumstats.gz"; mk_ovl "$L/ovl.gp2a.txt" pd_gp2 rls_akcimen "$S_AKC"
mk_info "$L/info.gp2d.txt" pd_gp2 85106 85106 "$M/pd_gp2_eur_n170212.sumstats.gz" rls_didriksen 10257 470725 "$M/rls_didriksen.sumstats.gz"; mk_ovl "$L/ovl.gp2d.txt" pd_gp2 rls_didriksen "$S_DID"
mk_info "$L/info.gp2s.txt" pd_gp2 85106 85106 "$M/pd_gp2_eur_n170212.sumstats.gz" rls 7248 19802 "$M/rls.sumstats.gz"; mk_ovl "$L/ovl.gp2s.txt" pd_gp2 rls "$S_SCH"
run gp2_tox3 g1000_eur "gp2_akc,$L/info.gp2a.txt,$L/ovl.gp2a.txt,pd_gp2,rls_akcimen;gp2_did,$L/info.gp2d.txt,$L/ovl.gp2d.txt,pd_gp2,rls_didriksen;gp2_sch,$L/info.gp2s.txt,$L/ovl.gp2s.txt,pd_gp2,rls" "$L/lava_pd_rls_gp2_tox3.csv" "$L/tox3.locfile"

echo "### 4. African ancestry, genome-wide ###"
mk_info "$L/info.afr.txt" pd_afr 1200 2445 "$M/pd_rizig_afr_n3645.sumstats.gz" rls_afr 2176 153313 "$M/rls_afr.sumstats.gz"; mk_ovl "$L/ovl.afr.txt" pd_afr rls_afr 0
run afr g1000_afr "afr,$L/info.afr.txt,$L/ovl.afr.txt,pd_afr,rls_afr" "$L/lava_pd_rls_afr.csv" "$LOC"
awk -F, 'END{print "tests:", NR-1}' "$L/lava_pd_rls_afr.csv"; awk -F, 'NR==1 || $0 ~ /52041337/' "$L/lava_pd_rls_afr.csv"
echo "### XANC_LAVA_TOX3_DONE ###"
