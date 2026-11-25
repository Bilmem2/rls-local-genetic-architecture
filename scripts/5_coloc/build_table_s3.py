#!/usr/bin/env python3
"""Supplementary Table S3: the 74 local genetic correlations carried to colocalization.
Inputs : results/known_vs_novel_blocks.csv   (blocks in hg19 and hg38; genome-wide-significance class)
         results/lava_test_family.tsv        (rho, P, Bonferroni P and BH q over the 2,739 tests)
         <analysis root>/results/coloc_screen_blocks.tsv  (coloc.abf over each LAVA block; coloc_screen_full.R)
Output : results/coloc_screen_full_supplementary.csv
Run from the repository root (as run_all.sh does)."""
import csv, io, os
ROOT = "/path/to/analysis"
RES = os.path.join(os.environ.get("REPO", "."), "results")
SCREEN = os.path.join(ROOT, "results", "coloc_screen_blocks.tsv")
code = {"anxiety": "anx", "bipolar": "bip", "chronotype": "chronotype", "daytime sleepiness": "daytime_sleepiness",
        "insomnia": "insomnia", "MDD": "mdd", "narcolepsy": "narcolepsy_ollila", "OSA": "osa", "PD": "pd_nalls",
        "PTSD": "ptsd", "RBD": "rbd", "schizophrenia": "scz", "long sleep": "sleepdur_eur_long",
        "short sleep": "sleepdur_eur_short"}
t1code = {"RLS": "rls_akcimen", "PD": "pd_nalls"}
gws = {"both": "both", "trait 1 only": "T1", "trait 2 only": "T2", "neither": "none", "neither (narcolepsy n/a)": "n/a"}
K = list(csv.DictReader(io.open(os.path.join(RES, "known_vs_novel_blocks.csv"), encoding="utf-8")))
U = {}
for r in csv.DictReader(io.open(os.path.join(RES, "lava_test_family.tsv"), encoding="utf-8"), delimiter="\t"):
    U[(r["phen1"], r["phen2"], r["chr"], r["start"], r["stop"])] = r
B = {r["k"]: r for r in csv.DictReader(io.open(SCREEN, encoding="utf-8"), delimiter="\t")}
assert len(K) == 74 and len(B) == 74
rows = []
for k in K:
    c, s, e = k["Region_hg19"].replace("chr", "").replace("-", ":").split(":")
    u = U[(t1code[k["Trait_1"]], code[k["Trait_2"]], c, s, e)]
    b = B["%s %s %s" % (k["Trait_1"], k["Trait_2"], k["Region_hg19"])]
    pp4, pp3 = float(b["PP4"]), float(b["PP3"])
    verdict = "COLOC" if pp4 > 0.90 else ("suggestive" if pp4 > 0.80 else "--")
    rows.append(dict(t1=k["Trait_1"], t2=k["Trait_2"], hg38=k["region_hg38"], hg19=k["Region_hg19"],
                     rho=float(u["rho"]), p=float(u["p"]), pb=float(u["p_bonf"]), q=float(u["q_BH"]),
                     nsnp=int(b["nSNP"]), pp3=pp3, pp4=pp4, lead=b["lead"], verdict=verdict, gws=gws[k["class"]],
                     pp4_strict=float(b["PP4_strict"]), pp4_lenient=float(b["PP4_lenient"])))
assert all(r["q"] < 0.05 for r in rows) and sum(r["pb"] < 0.05 for r in rows) == 14
rows.sort(key=lambda r: (-round(r["pp4"], 3), -round(r["pp3"], 3)))
with io.open(os.path.join(RES, "coloc_screen_full_supplementary.csv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Trait_1", "Trait_2", "Region_hg19", "region_hg38", "LAVA_rho", "LAVA_P", "p_bonferroni", "q_BH",
                "survives_bonferroni", "nSNP", "PP.H3", "PP.H4", "PP.H4_p12_1e-6", "PP.H4_p12_5e-5", "lead_SNP", "verdict", "GWS"])
    for r in rows:
        w.writerow([r["t1"], r["t2"], r["hg19"], r["hg38"], "%.3f" % r["rho"], "%.3g" % r["p"], "%.3g" % r["pb"],
                    "%.3g" % r["q"], "yes" if r["pb"] < 0.05 else "no", r["nsnp"], "%.3f" % r["pp3"], "%.3f" % r["pp4"],
                    "%.3f" % r["pp4_strict"], "%.3f" % r["pp4_lenient"], r["lead"], r["verdict"], r["gws"]])
print("rows", len(rows), "| COLOC", sum(r["verdict"] == "COLOC" for r in rows),
      "| suggestive", sum(r["verdict"] == "suggestive" for r in rows))
