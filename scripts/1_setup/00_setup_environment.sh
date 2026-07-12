#!/usr/bin/env bash
# Reproducible record of the environment build (Phase 0).
# Not meant to be re-run blindly; documents exactly how the toolchain was created.
set -euo pipefail
source "$(conda info --base 2>/dev/null || echo "$HOME/miniconda3")/etc/profile.d/conda.sh"

# 1. Dedicated, version-pinned LDSC env (existing envs carry numpy 2.x / pandas 3.x,
#    which break LDSC v1.0.1). numpy 1.23.5 keeps np.int alive; pandas 1.4.4 is the
#    last branch the .ix->.loc patches target cleanly.
mamba create -y -n ldsc -c conda-forge -c bioconda \
  python=3.9 numpy=1.23.5 scipy=1.9.3 pandas=1.4.4 bitarray=2.6.0 \
  plink plink2 git
# NOTE: /usr/bin/plink on this box is PuTTY's plink (v0.81), NOT genetics PLINK.
# bioconda 'plink' (v1.90b7.7) + 'plink2' installed into the env instead.

# 2. LDSC v1.0.1, python-3 patched.
conda activate ldsc
cd /path/to/analysis/tools
git clone --depth 1 https://github.com/bulik/ldsc.git
cd ldsc
python -m lib2to3 -w -n ldsc.py munge_sumstats.py ldscore/*.py        # py2 -> py3
for f in ldsc.py munge_sumstats.py ldscore/*.py; do                   # fix mixed tabs/spaces
  expand -i -t 8 "$f" > "$f.tmp" && mv "$f.tmp" "$f"; done
# pandas>=1.0 API fixes:
sed -i "s/set_option('precision', 4)/set_option('display.precision', 4)/" ldsc.py
sed -i 's/ref_ld\.ix\[:, 1:\]\.var()/ref_ld.iloc[:, 1:].var()/' ldscore/sumstats.py
sed -i 's/ref_ld = ref_ld\.ix\[:, ii_snp\]/ref_ld = ref_ld.loc[:, ii_snp]/' ldscore/sumstats.py
sed -i 's/sumstats = sumstats\.ix\[ii, :\]/sumstats = sumstats.loc[ii, :]/g' ldscore/sumstats.py
sed -i 's/\.ix\[:,1:\])/.iloc[:,1:])/' ldscore/sumstats.py
sed -i 's/sumstats\.as_matrix(columns=ref_ld_cnames)/sumstats[ref_ld_cnames].values/' ldscore/sumstats.py
sed -i 's/df = df\.ix\[df\.SNP\.isin(print_snps\.SNP),:\]/df = df.loc[df.SNP.isin(print_snps.SNP),:]/' ldsc.py
sed -i 's/t = df\.ix\[:,4:\]\.describe()/t = df.iloc[:,4:].describe()/' ldsc.py
sed -i 's/t\.ix\[1:,:\]/t.iloc[1:,:]/' ldsc.py
sed -i 's/df\.ix\[:,4:\]\.corr()/df.iloc[:,4:].corr()/' ldsc.py
sed -i 's/np\.linalg\.cond(df\.ix\[:,5:\])/np.linalg.cond(df.iloc[:,5:])/' ldsc.py
# py3 I/O fixes surfaced by the end-to-end smoke test:
sed -i "s/self.log_fh = open(fh, 'wb')/self.log_fh = open(fh, 'w')/" ldsc.py   # log handle was binary
sed -i "s/'.M','wb')/'.M','w')/;s/'.M_5_50','wb')/'.M_5_50','w')/" ldsc.py       # .M writers binary
sed -i 's/traceback\.format_exc(ex)/traceback.format_exc()/' munge_sumstats.py   # format_exc takes no exc arg in py3
# munge reads gzip header in binary by default -> str/bytes TypeError in read_header; open compressed in text mode:
sed -i 's/openfunc = gzip.open/openfunc = lambda f: gzip.open(f, "rt")/' munge_sumstats.py
sed -i 's/openfunc = bz2.BZ2File/openfunc = lambda f: bz2.open(f, "rt")/' munge_sumstats.py

# Validated end-to-end: munge x2 (1.217M HM3 SNPs each) + rg ran clean on synthetic
# data (rg=nan / "h2 out of bounds" is the CORRECT result for signal-free random Z).

# 3. LAVA 0.1.5 was already installed in system R 4.3.3 (user library); confirmed loading.
echo "environment build recorded."
