#!/usr/bin/env python3
# Figure 5 - TOX3/CASC16 (chr16q12, LAVA block 2134) local rg forest: RLS x PD across PD GWAS
# and ancestries.
#
# Every estimate is read from results/tox3_crossancestry_local_rg.csv (shipped with this
# repository), so the figure cannot drift away from the reported numbers. Rows carrying a
# confidence interval are gated LAVA estimates; rows without one are ungated point estimates
# and are drawn as open diamonds. This script fixes only labels, colours and order.
import csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = "/path/to/analysis"
CSV  = ROOT + "/results/tox3_crossancestry_local_rg.csv"
OUT  = ROOT + "/results/figures/fig_forest_TOX3.png"

# display names; unknown keys raise rather than silently mislabelling a dataset
PD_SHORT  = {"Loesch 2021": "Loesch", "Rizig 2023": "Rizig", "FinnGen R11": "FinnGen",
             "GP2 (Leonard 2025)": "GP2", "Nalls 2019": "Nalls"}
RLS_SHORT = {"RLS-AMR (Akcimen)": "RLS-AMR", "RLS-AFR (Akcimen)": "RLS-AFR",
             "Schormair 2024": "Schormair", "Didriksen 2020": "Didriksen", "Akcimen 2026": "Akçimen"}

D = []
for r in csv.DictReader(open(CSV, encoding="utf-8")):
    pd_s, rls_s, anc = PD_SHORT[r["pd_gwas"]], RLS_SHORT[r["rls_gwas"]], r["ancestry"]
    ci = bool(r["ci_low"].strip())
    D.append((f"{anc}: {pd_s} × {rls_s}", float(r["local_rg_rho"]), float(r["ci_low"]) if ci else None,
              float(r["ci_high"]) if ci else None, anc if anc != "EUR" else f"EUR-{pd_s}", ci))
order = ["EUR: Nalls × Akçimen", "EUR: Nalls × Didriksen", "EUR: Nalls × Schormair", "EUR: GP2 × Akçimen",
         "EUR: GP2 × Didriksen", "EUR: GP2 × Schormair", "EUR: FinnGen × Akçimen", "EUR: FinnGen × Didriksen",
         "EUR: FinnGen × Schormair", "AFR: Rizig × RLS-AFR", "AMR: Loesch × RLS-AMR"]
D.sort(key=lambda d: order.index(d[0])); assert len(D) == 11 and all(d[1] < 0 for d in D)
col = {"EUR-Nalls": "#1f4e79", "EUR-GP2": "#2a8a8a", "EUR-FinnGen": "#7aa6c2", "AFR": "#b1322e", "AMR": "#c98a2e"}
fig, ax = plt.subplots(figsize=(8.4, 5.6))
for i, (lab, rho, lo, hi, g, ci) in enumerate(D):
    c = col[g]
    if ci:
        ax.plot([lo, hi], [i, i], color=c, lw=2.2, solid_capstyle="round")
        ax.scatter([rho], [i], s=70, color=c, zorder=3, edgecolor="white", linewidth=0.8, marker="o")
    else:
        ax.scatter([rho], [i], s=46, facecolor="white", edgecolor=c, zorder=3, linewidth=1.6, marker="D")
ax.axvline(0, color="k", lw=1, ls="--")
ax.set_yticks(range(len(D))); ax.set_yticklabels([d[0] for d in D], fontsize=9); ax.invert_yaxis()
ax.set_xlabel("Local genetic correlation $\\rho$  (RLS × PD at TOX3/CASC16)")
ax.set_xlim(-1.05, 0.35)
ax.set_title("TOX3/CASC16 (chr16q12): the negative RLS–PD local rg is directionally\n"
             f"consistent across PD GWAS and ancestries (all {len(D)} estimates negative)", fontsize=11)
ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor="#1f4e79", markersize=8, label="filled = gated LAVA estimate, 95% CI"),
                   Line2D([0], [0], marker="D", color="w", markerfacecolor="white", markeredgecolor="#7aa6c2", markeredgewidth=1.6, markersize=7, label="open diamond = ungated point estimate")],
          loc="lower left", fontsize=8, frameon=True)
ax.grid(axis="x", alpha=0.25)
plt.tight_layout(); plt.savefig(OUT, dpi=300, bbox_inches="tight"); plt.close()
print("[saved]", OUT)
