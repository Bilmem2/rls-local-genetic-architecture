#!/usr/bin/env python3
# TOX3/CASC16 (chr16q12, block 2134) local rg forest: RLS x PD across datasets + ancestries.
# CIs from gated LAVA CSVs (Nalls/GP2/AFR); FinnGen + AMR shown as point estimates (from the 11-pair triage).
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
# label, rho, lo, hi, group, has_ci
D=[
 ("AMR: Loesch × RLS-AMR",      -0.333, None, None, "AMR", False),
 ("AFR: Rizig × RLS-AFR",       -0.734, -1.000, -0.419, "AFR", True),
 ("EUR: FinnGen × Schormair",   -0.025, None, None, "EUR-FinnGen", False),
 ("EUR: FinnGen × Didriksen",   -0.329, None, None, "EUR-FinnGen", False),
 ("EUR: FinnGen × Akçimen",-0.098, None, None, "EUR-FinnGen", False),
 ("EUR: GP2 × Schormair",       -0.306, -0.496, -0.122, "EUR-GP2", True),
 ("EUR: GP2 × Didriksen",       -0.556, -0.733, -0.378, "EUR-GP2", True),
 ("EUR: GP2 × Akçimen",    -0.490, -0.649, -0.333, "EUR-GP2", True),
 ("EUR: Nalls × Schormair",     -0.452, -0.762, -0.184, "EUR-Nalls", True),
 ("EUR: Nalls × Didriksen",     -0.662, -0.998, -0.410, "EUR-Nalls", True),
 ("EUR: Nalls × Akçimen",  -0.530, -0.826, -0.292, "EUR-Nalls", True),
]
col={"EUR-Nalls":"#1f4e79","EUR-GP2":"#2a8a8a","EUR-FinnGen":"#7aa6c2","AFR":"#b1322e","AMR":"#c98a2e"}
fig,ax=plt.subplots(figsize=(8.4,5.6))
for i,(lab,rho,lo,hi,g,ci) in enumerate(D):
    c=col[g]
    if ci: ax.plot([lo,hi],[i,i],color=c,lw=2.2,solid_capstyle="round")
    ax.scatter([rho],[i],s=70 if ci else 46,color=c,zorder=3,
               edgecolor="white",linewidth=0.8,marker="o" if ci else "D")
ax.axvline(0,color="k",lw=1,ls="--")
ax.set_yticks(range(len(D))); ax.set_yticklabels([d[0] for d in D],fontsize=9)
ax.set_xlabel("Local genetic correlation $\\rho$  (RLS × PD at TOX3/CASC16)")
ax.set_xlim(-1.05,0.35)
ax.set_title("TOX3/CASC16 (chr16q12): RLS–PD negative local rg replicates\nacross PD GWAS and ancestries (all 11 estimates negative)",fontsize=11)
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],marker="o",color="w",markerfacecolor="#1f4e79",markersize=8,label="filled = 95% CI (gated LAVA)"),
     Line2D([0],[0],marker="D",color="w",markerfacecolor="#7aa6c2",markersize=7,label="open ◆ = point est. (triage, underpowered)")]
ax.legend(handles=leg,loc="lower left",fontsize=8,frameon=True)
ax.grid(axis="x",alpha=0.25)
plt.tight_layout(); plt.savefig("/path/to/analysis/results/lava/fig_forest_TOX3.png",dpi=300,bbox_inches="tight")
print("[saved] fig_forest_TOX3.png")
