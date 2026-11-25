#!/usr/bin/env bash
# Effective sample size from standard errors, for a GWAS file without a per-variant N:
#   k = median over common variants (0.05 < f < 0.95) of 1 / (2 f (1 - f) SE^2); for a case-control GWAS
#   on the log-odds scale, Neff = 4 / (1/cases + 1/controls) is approximately 4k.
# Checked on FinnGen R11 PD (G6_PARKINSON, 5,150 cases / 448,583 controls: Neff 20,365), then applied to the
# GP2 European meta-analysis (used in xanc_lava_tox3.sh). Every fifth line is read.
ROOT=/path/to/analysis; G=$ROOT/data/gwas
est(){ # label, reader command, SE column, effect-allele-frequency column
  local tag=$1 cmd=$2 sc=$3 fc=$4
  eval "$cmd" 2>/dev/null | awk -v sc=$sc -v fc=$fc -v tag="$tag" 'NR>1{ f=$fc+0; s=$sc+0; if(f>0.05&&f<0.95&&s>0){ v[++n]=1/(2*f*(1-f)*s*s) } }
    END{ asort(v); printf "%-30s SNPs=%d  k=%.0f  4k=%.0f\n", tag, n, v[int((n+1)/2)], 4*v[int((n+1)/2)] }'
}
est "FinnGen R11 PD (Neff 20,365)" "zcat $G/pd/finngen_R11_PD.gz | awk 'NR%5==1'" 10 11
est "GP2 European meta-analysis" "unzip -p $G/pd/GP2_euro_ancestry_meta_analysis_2024.zip GP2_euro_ancestry_meta_analysis_2024/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt.gz | zcat | awk 'NR%5==1'" 10 6
