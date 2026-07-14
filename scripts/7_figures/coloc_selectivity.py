#!/usr/bin/env python3
# Figure 3 - colocalization filters local genetic correlation: every significant local rg in the
# panel plotted against its colocalization posterior, showing how few share a causal variant.
#
# Points and the counts quoted in the title are read from results/coloc_screen_full_supplementary.csv
# (shipped with this repository), so the figure cannot drift away from the reported numbers. Only
# the label text and colours are fixed here.
#
# Requires: matplotlib, adjustText  (pip install matplotlib adjustText)
import csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

ROOT = "/path/to/analysis"
CSV  = ROOT + "/results/coloc_screen_full_supplementary.csv"
OUT  = ROOT + "/results/figures/coloc_selectivity.png"

rows = list(csv.DictReader(open(CSV, encoding="utf-8")))

# (region start, trait keyword) -> italic label for the loci worth naming on the plot
LABELS = [("2:65938003",  "insomnia",    "$MEIS1$ $\\times$ insomnia"),
          ("2:65938003",  "short sleep", "$MEIS1$ $\\times$ short-sleep"),
          ("2:65938003",  "chronotype",  "$MEIS1$ $\\times$ chronotype"),
          ("9:8262304",   "short sleep", "$PTPRD$ $\\times$ short-sleep"),
          ("15:47693317", "MDD",         "$SEMA6D$ $\\times$ MDD"),
          ("16:52041337", "PD",          "$TOX3$ $\\times$ PD"),
          ("18:44300146", "daytime",     "$SKOR2$ $\\times$ daytime"),
          ("1:154685546", "RBD",         "$GBA$ reg. $\\times$ RBD")]

def labfor(r):
    reg = r["Region_hg19"].replace("chr", "")
    for start, tr, lab in LABELS:
        if reg.startswith(start) and tr.split()[0].lower() in r["Trait"].lower():
            return lab
    return None

def style(v):
    if v == "COLOC":      return "#c0392b", 105, 4
    if v == "suggestive": return "#e08e0b", 82, 4
    return "#b3b3b3", 26, 2

n_total = len(rows)
n_coloc = sum(1 for r in rows if r["verdict"] == "COLOC")
n_genes = len({lab.split("$")[1] for lab in (labfor(r) for r in rows if r["verdict"] == "COLOC") if lab})

fig, ax = plt.subplots(figsize=(9.2, 6.2))
# diffuse points first (behind), then labelled points on top with a dark edge
for r in sorted(rows, key=lambda r: r["verdict"] != "no"):
    c, s, z = style(r["verdict"])
    lab = labfor(r)
    ax.scatter(float(r["LAVA_rho"]), float(r["PP.H4"]), s=s, c=c, marker="o",
               edgecolor=("#333333" if lab else "white"), linewidth=(0.9 if lab else 0.5), zorder=z)
ax.axhline(0.9, ls="--", color="#c0392b", lw=1, alpha=0.75)
ax.text(-1.05, 0.925, "PP4 = 0.9 (coloc threshold)", color="#c0392b", fontsize=8)

texts = [ax.text(float(r["LAVA_rho"]), float(r["PP.H4"]), labfor(r), fontsize=9, zorder=6)
         for r in rows if labfor(r)]
adjust_text(texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="0.5", lw=0.7),
            expand=(1.6, 2.0), force_text=(0.6, 1.0), force_static=(0.3, 0.5),
            only_move={"text": "xy", "static": "xy"})

ax.set_xlabel("LAVA local genetic correlation ($\\rho$)", fontsize=11)
ax.set_ylabel("Colocalization posterior PP4 (shared causal variant)", fontsize=11)
ax.set_title("Colocalization filters local genetic correlation:\n"
             f"{n_coloc} of {n_total} panel hits share a causal variant  "
             "($MEIS1$ shared by 3 sleep traits)", fontsize=11.5)
leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#b3b3b3", markersize=7,
              label="diffuse (no shared causal variant)"),
       Line2D([0], [0], marker="o", color="w", markerfacecolor="#e08e0b", markersize=8,
              label="suggestive (0.8 < PP4 < 0.9)"),
       Line2D([0], [0], marker="o", color="w", markerfacecolor="#c0392b", markersize=9,
              label="colocalizing (PP4 > 0.9)")]
ax.legend(handles=leg, loc="lower center", fontsize=8.5, frameon=True, ncol=1)
ax.set_xlim(-1.12, 1.12); ax.set_ylim(-0.05, 1.08); ax.grid(axis="y", alpha=0.2)
plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches="tight")
print(f"[saved] {OUT}  ({n_coloc} colocalizing of {n_total} hits, at {n_genes} genes)")
