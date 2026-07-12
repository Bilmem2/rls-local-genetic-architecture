#!/usr/bin/env python3
# Annotate rsIDs onto hg38 PD sumstats that lack them (GP2 EUR, Loesch AMR), using the HM3 hg38 map.
# Match on (chr, pos_hg38) + allele-set; output clean TSV (SNP A1 A2 BETA P) for LDSC munge.
# N: Loesch has per-SNP n.obs (median used); GP2 has none -> effective N derived as median 1/(2f(1-f)SE^2).
import gzip, subprocess, statistics
ROOT="/path/to/analysis"
MAP=f"{ROOT}/data/reference/xanc_map/hm3_hg38_map.tsv"
PD=f"{ROOT}/data/gwas/pd"; OUT=f"{ROOT}/results/munged"

m={}
with open(MAP) as f:
    for line in f:
        c,p,rs,a1,a2=line.rstrip("\n").split("\t")
        m[(c,p)]=(rs,a1.upper(),a2.upper())
print(f"[map] {len(m)} HM3 hg38 entries", flush=True)

def annotate(name, it, col, nmode, outpath):
    nvals=[]; kept=0; tot=0
    with open(outpath,"w") as o:
        o.write("SNP\tA1\tA2\tBETA\tP\n")
        first=True
        for line in it:
            if first: first=False; continue
            f=line.rstrip("\n").split("\t")
            try:
                c=f[col["chr"]].replace("chr",""); p=f[col["pos"]]
                ea=f[col["ea"]].upper(); oa=f[col["oa"]].upper()
                beta=f[col["beta"]]; pv=f[col["p"]]
            except Exception:
                continue
            tot+=1
            hit=m.get((c,p))
            if not hit: continue
            rs,A1,A2=hit
            if {ea,oa}!={A1,A2}: continue
            o.write(f"{rs}\t{ea}\t{oa}\t{beta}\t{pv}\n"); kept+=1
            if nmode=="nobs":
                try: nvals.append(float(f[col["nobs"]]))
                except Exception: pass
            else:
                try:
                    se=float(f[col["se"]]); fr=float(f[col["eaf"]])
                    if 0.01<fr<0.99 and se>0: nvals.append(1.0/(2*fr*(1-fr)*se*se))
                except Exception: pass
    Nuse=int(statistics.median(nvals)) if nvals else 0
    print(f"[{name}] scanned={tot} matched(rsID+alleles)={kept} N({nmode} median)={Nuse}", flush=True)
    return Nuse

# GP2 (inside zip): chromosome[0] base_pair_position[1] SNP_ID[2] effect_allele[3] other_allele[4] eaf[5] n_datasets[6] p_value[7] beta[8] SE[9]
gp2_zip=f"{PD}/GP2_euro_ancestry_meta_analysis_2024.zip"
gp2_inner="GP2_euro_ancestry_meta_analysis_2024/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt.gz"
pr=subprocess.Popen(["bash","-c",f"unzip -p '{gp2_zip}' '{gp2_inner}' | zcat"], stdout=subprocess.PIPE, text=True)
gp2N=annotate("GP2", pr.stdout, {"chr":0,"pos":1,"ea":3,"oa":4,"eaf":5,"p":7,"beta":8,"se":9}, "derive", f"{OUT}/pd_gp2_eur.clean.tsv")
pr.wait()

# Loesch: #SNP[0] variant.id[1] chr[2] pos[3] ref[4] alt[5] allele.index[6] n.obs[7] freq[8] ... Score.pval[14] beta[15]
lo=f"{PD}/release5_11052023_summary_statistics_Loesch_et_al_2021_Latam_no23andme_hg38.txt.gz"
loN=annotate("Loesch", gzip.open(lo,"rt"), {"chr":2,"pos":3,"oa":4,"ea":5,"nobs":7,"p":14,"beta":15}, "nobs", f"{OUT}/pd_loesch_amr.clean.tsv")

with open(f"{OUT}/pd_xanc.N","w") as f:
    f.write(f"GP2 {gp2N}\nLOESCH {loN}\n")
print("[done] clean TSVs + pd_xanc.N written", flush=True)
