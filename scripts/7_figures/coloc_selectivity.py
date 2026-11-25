#!/usr/bin/env python3
# Figure 3 - local genetic correlation against colocalization posterior for the 74 local correlations
# carried to colocalization.
#
# Points and the counts in the title are read from results/coloc_screen_full_supplementary.csv
# (shipped with this repository; written by scripts/5_coloc/build_table_s3.py). Only label text,
# label positions and colours are fixed here.
import csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = "/path/to/analysis"
CSV  = ROOT + "/results/coloc_screen_full_supplementary.csv"
OUT  = ROOT + "/results/figures/coloc_selectivity.png"

rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
assert len(rows) == 74
LAB = {("RLS", "chronotype", "2:65710869-66996895"): ("$MEIS1$ × chronotype", (-0.72, 1.045)),
       ("RLS", "PD", "16:52007425-53359970"): ("$TOX3$ × PD", (-0.43, 0.965)),
       ("RLS", "short sleep", "2:65710869-66996895"): ("$MEIS1$ × short sleep", (0.06, 1.045)),
       ("RLS", "OSA", "2:158184468-159343868"): ("$CCDC148/PKP4$ × OSA", (0.02, 0.855)),
       ("RLS", "daytime sleepiness", "18:46720183-48192121"): ("$SKOR2$ × daytime\nsleepiness", (0.755, 0.855)),
       ("RLS", "insomnia", "2:65710869-66996895"): ("$MEIS1$ × insomnia", (0.66, 1.045)),
       ("RLS", "MDD", "15:47401120-48967809"): ("$SEMA6D$ × MDD", (0.79, 0.935)),
       ("PD", "RBD", "1:154713070-156844053"): ("$GBA1$ × RBD", (0.68, 0.755))}
fig, ax = plt.subplots(figsize=(10, 6.7))
n = {"COLOC": 0, "suggestive": 0}
for r in rows:
    x, y, v = float(r["LAVA_rho"]), float(r["PP.H4"]), r["verdict"]
    key = (r["Trait_1"], r["Trait_2"], r["region_hg38"])
    if v == "--":
        assert key not in LAB
        ax.scatter(x, y, s=34, color="#b3b3b3", edgecolor="white", linewidth=0.6, zorder=2)
        continue
    n[v] += 1
    col = "#c0392b" if v == "COLOC" else "#e08e0b"
    ax.scatter(x, y, s=150, color=col, edgecolor="#333333", linewidth=1.0, zorder=4)
    lab, (tx, ty) = LAB[key]
    ax.annotate(lab, (x, y), (tx, ty), fontsize=11, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color="#777777", lw=0.9, shrinkA=2, shrinkB=7))
assert n == {"COLOC": 7, "suggestive": 1}, n
ax.axhline(0.9, color="#c0392b", ls="--", lw=1.1, alpha=0.8)
ax.text(-1.05, 0.925, "PP4 = 0.9 (coloc threshold)", color="#c0392b", fontsize=9)
ax.set_xlim(-1.12, 1.12); ax.set_ylim(-0.05, 1.08)
ax.set_xlabel(r"LAVA local genetic correlation ($\rho$)", fontsize=12)
ax.set_ylabel("Colocalization posterior PP4 (shared causal variant)", fontsize=12)
ax.set_title("Colocalization of the 74 screened local genetic correlations:\n"
             "7 colocalize (PP4 > 0.9), three of them at $MEIS1$", fontsize=13)
ax.grid(axis="y", alpha=0.25)
ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor="#b3b3b3", markersize=8, label="not colocalizing (PP4 < 0.8)"),
                   Line2D([0], [0], marker="o", color="w", markerfacecolor="#e08e0b", markersize=10, label="suggestive (0.8 < PP4 < 0.9)"),
                   Line2D([0], [0], marker="o", color="w", markerfacecolor="#c0392b", markersize=10, label="colocalizing (PP4 > 0.9)")],
          loc="lower center", fontsize=10, frameon=True)
plt.tight_layout(); plt.savefig(OUT, dpi=270); plt.close()
print("[saved]", OUT)
