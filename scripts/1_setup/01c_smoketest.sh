#!/usr/bin/env bash
# End-to-end py3 LDSC smoke test: synthetic munge x2 + rg, using Phase3 reference.
set -uo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"
conda activate ldsc
ROOT=/path/to/analysis
REF=$ROOT/data/ldsc_ref
LDSC=$ROOT/tools/ldsc
REFLD="$REF/LDscore/LDscore."
WLD="$REF/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
SMK=$ROOT/results/_smoketest
mkdir -p "$SMK"; cd "$SMK"

python - <<'PY'
import numpy as np, pandas as pd
from scipy.stats import norm
snp = pd.read_csv('/path/to/analysis/data/ldsc_ref/w_hm3.snplist', sep=r'\s+')
snp = snp.sample(n=120000, random_state=1).reset_index(drop=True)
rng = np.random.default_rng(0)
for tag in ('a','b'):
    z = rng.standard_normal(len(snp))
    p = 2*norm.sf(np.abs(z))
    pd.DataFrame({'SNP':snp.SNP,'A1':snp.A1,'A2':snp.A2,'N':100000,
                  'Z':z,'P':p}).to_csv(f'synth_{tag}.txt',sep='\t',index=False)
print('wrote synthetic sumstats rows=',len(snp))
PY

for tag in a b; do
  python "$LDSC/munge_sumstats.py" --sumstats synth_$tag.txt --N 100000 \
    --merge-alleles "$REF/w_hm3.snplist" --signed-sumstats Z,0 \
    --snp SNP --a1 A1 --a2 A2 --p P --out synth_$tag > munge_$tag.log 2>&1 \
    && echo "munge_$tag OK ($(zcat synth_$tag.sumstats.gz | tail -n +2 | wc -l) SNPs)" \
    || { echo "munge_$tag FAIL"; tail -20 munge_$tag.log; }
done

if [ -f synth_a.sumstats.gz ] && [ -f synth_b.sumstats.gz ]; then
  python "$LDSC/ldsc.py" --rg synth_a.sumstats.gz,synth_b.sumstats.gz \
    --ref-ld-chr "$REFLD" --w-ld-chr "$WLD" --out smoke_rg > smoke_rg.log 2>&1 \
    && { echo "=== rg OK ==="; grep -E "Total Observed scale h2|Genetic Correlation:|Z-score:|P:|Intercept:" smoke_rg.log; } \
    || { echo "=== rg FAIL ==="; tail -30 smoke_rg.log; }
fi
echo "=== SMOKE_DONE ==="
