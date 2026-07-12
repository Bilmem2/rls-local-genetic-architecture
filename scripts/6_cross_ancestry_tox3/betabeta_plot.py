#!/usr/bin/env python3
# Beta/beta plots at TOX3/CASC16 (Fulya's suggestion): Nalls PD vs each RLS GWAS.
# rsID-matched, alleles harmonised to PD effect allele (RLS beta flipped if swapped).
# Plotted on locus-ASSOCIATED SNPs (signal in either trait) like Fulya's figure; all-SNP r also reported.
# This is a DESCRIPTIVE visualisation; the rigorous estimate is the LAVA local rho.
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

pdf=pd.read_csv("pd_nalls_tox3.txt",sep=r"\s+",header=None,names=["rsid","ea","oa","beta_pd","p_pd"])
for c in ["beta_pd","p_pd"]: pdf[c]=pd.to_numeric(pdf[c],errors="coerce")
pdf=pdf.dropna(subset=["rsid","ea","oa","beta_pd","p_pd"]).drop_duplicates("rsid")
pdf["ea"]=pdf.ea.astype(str).str.upper(); pdf["oa"]=pdf.oa.astype(str).str.upper()

rls={"Akçimen":("rls_akcimen_chr16.txt","beta",-0.53),
     "Didriksen":("rls_didriksen_chr16.txt","or",-0.66),
     "Schormair":("rls_schormair_chr16.txt","beta",-0.45)}

fig,axes=plt.subplots(1,3,figsize=(15,4.8))
for ax,(name,(f,kind,rho)) in zip(axes,rls.items()):
    r=pd.read_csv(f,sep=r"\s+",header=None,names=["rsid","ea","oa","val","p_rls"])
    for c in ["val","p_rls"]: r[c]=pd.to_numeric(r[c],errors="coerce")
    r=r.dropna(subset=["rsid","ea","oa","val","p_rls"]).drop_duplicates("rsid")
    r["ea"]=r.ea.astype(str).str.upper(); r["oa"]=r.oa.astype(str).str.upper()
    if kind=="or":
        r=r[r.val>0]; r["beta_rls"]=np.log(r.val)
    else:
        r["beta_rls"]=r.val
    m=pdf.merge(r,on="rsid",suffixes=("_pd","_rls"))
    same=(m.ea_pd==m.ea_rls)&(m.oa_pd==m.oa_rls)
    swap=(m.ea_pd==m.oa_rls)&(m.oa_pd==m.ea_rls)
    m=m[same|swap].copy(); m.loc[swap,"beta_rls"]=-m.loc[swap,"beta_rls"]
    m=m[np.isfinite(m.beta_rls)&np.isfinite(m.beta_pd)&np.isfinite(m.p_pd)&np.isfinite(m.p_rls)]
    r_all=np.corrcoef(m.beta_rls,m.beta_pd)[0,1] if len(m)>2 else float("nan")
    # locus-associated SNPs (signal in either trait)
    minp=np.minimum(m.p_pd.clip(lower=1e-300),m.p_rls.clip(lower=1e-300))
    thr=1e-4; a=m[minp<thr]
    if len(a)<25: thr=1e-3; a=m[minp<thr]
    x=a.beta_rls.values; y=a.beta_pd.values
    rr=np.corrcoef(x,y)[0,1] if len(a)>2 else float("nan")
    w=-np.log10(np.minimum(a.p_pd.clip(lower=1e-300),a.p_rls.clip(lower=1e-300)).values) if len(a) else np.array([1.0])
    s=8+26*(w-w.min())/(w.max()-w.min()+1e-9)
    ax.axhline(0,ls="--",c="k",lw=.7); ax.axvline(0,ls="--",c="k",lw=.7)
    ax.scatter(x,y,s=s,c="#2a8a8a",alpha=.72,edgecolor="none")
    if len(a)>2:
        b1,b0=np.polyfit(x,y,1); xs=np.array([x.min(),x.max()]); ax.plot(xs,b0+b1*xs,c="#b1322e",lw=2)
    ax.set_title(f"TOX3/CASC16 — Nalls × {name}",fontsize=11,style="italic")
    ax.set_xlabel("BETA in RLS"); ax.set_ylabel("BETA in PD (Nalls 2019)")
    ax.text(.96,.96,f"beta–beta r = {rr:.2f}  (assoc, n={len(a)})\nall-SNP r = {r_all:.2f}  (n={len(m)})\nLAVA ρ = {rho:.2f}",
            transform=ax.transAxes,ha="right",va="top",fontsize=8.5,fontweight="bold")
    print(f"{name:<10} region n={len(m):>4} | assoc(minp<{thr:g}) n={len(a):>3} | beta-beta r(assoc)={rr:.3f} | r(all)={r_all:.3f} | LAVA rho={rho}")
plt.suptitle("TOX3/CASC16 (chr16q12) — per-SNP effect concordance, RLS vs PD (negative = opposite directions)",fontsize=12,y=1.02)
plt.tight_layout(); plt.savefig("fig_betabeta_TOX3.png",dpi=150,bbox_inches="tight")
print("[saved] fig_betabeta_TOX3.png")
