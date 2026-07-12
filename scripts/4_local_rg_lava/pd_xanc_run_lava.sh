#!/usr/bin/env bash
# =============================================================================
# FINAL cross-ancestry + newer-EUR-PD LAVA (genome-wide). RUN IN YOUR OWN WSL TERMINAL (~1-2 h total).
#   A3  GP2-2024 EUR PD x {Akçimen,Didriksen,Schormair} RLS  (g1000_eur) — newer/independent EUR PD, beyond Nalls/TOX3
#   AFR Rizig-2023 PD   x RLS-AFR                            (g1000_afr)
#   AMR Loesch-2021 PD  x RLS-AMR                            (g1000_amr) — very low power (N~1.5k), expect null
#   Each = pd_05_lava.R genome-wide scan (2495 blocks; both-trait univ gate p<0.05/2495; Bonferroni+BH).
#   Identity sample-overlap (independent-cohort assumption); EUR-hg19 blocks (rsID-matched, build-agnostic).
#   Smoke-tested @ TOX3:  GP2 rho=-0.49 (p2e-8) | AFR rho=-0.73 (p3e-5) | AMR rho=-0.33 (p0.20, underpowered).
# =============================================================================
set -uo pipefail
ROOT=/path/to/analysis
L="$ROOT/results/lava"; M="$ROOT/results/munged"; mkdir -p "$L" "$ROOT/results/_logs"
LOC="$ROOT/data/reference/lava_blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
RS=/usr/bin/Rscript
FILT='Extracting|Aligning|Negative variance|out of bounds|Removing|Reading|Processing input|treated as BINARY|^\|'
mk_info(){ printf 'phenotype cases controls filename\n%s %s %s %s\n%s %s %s %s\n' "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" > "$1"; }
mk_ovl(){ printf '%s %s\n%s 1 0\n%s 0 1\n' "$2" "$3" "$2" "$3" > "$1"; }
runscan(){ # 1 tag  2 ref-basename  3 LAVA_PAIRS  4 out.csv
  export LAVA_REF="$ROOT/data/reference/$2"; export LAVA_LOCF="$LOC"; export LAVA_OUT="$4"; export LAVA_PAIRS="$3"
  echo "### $1  (ref=$2)  -> $4 ###"; cd "$ROOT"
  stdbuf -oL -eL "$RS" scripts/pd_05_lava.R 2>&1 | tee "$ROOT/results/_logs/lava_$1.log" | grep --line-buffered -vE "$FILT"
}
echo "### CROSS-ANCESTRY + GP2 LAVA  $(date -u +%FT%TZ) ###"

# ---- A3: GP2 EUR PD x 3 RLS (g1000_eur) ----
mk_info "$L/info.gp2_akc.txt" pd_gp2 21503 21503 "$M/pd_gp2_eur.sumstats.gz" rls_akcimen   36993 639182 "$M/rls_akcimen.sumstats.gz";   mk_ovl "$L/ovl.gp2_akc.txt" pd_gp2 rls_akcimen
mk_info "$L/info.gp2_did.txt" pd_gp2 21503 21503 "$M/pd_gp2_eur.sumstats.gz" rls_didriksen 10257 470725 "$M/rls_didriksen.sumstats.gz"; mk_ovl "$L/ovl.gp2_did.txt" pd_gp2 rls_didriksen
mk_info "$L/info.gp2_sch.txt" pd_gp2 21503 21503 "$M/pd_gp2_eur.sumstats.gz" rls           7248  19802  "$M/rls.sumstats.gz";           mk_ovl "$L/ovl.gp2_sch.txt" pd_gp2 rls
runscan A3_gp2 g1000_eur "gp2_akcimen,$L/info.gp2_akc.txt,$L/ovl.gp2_akc.txt,pd_gp2,rls_akcimen;gp2_didriksen,$L/info.gp2_did.txt,$L/ovl.gp2_did.txt,pd_gp2,rls_didriksen;gp2_schormair,$L/info.gp2_sch.txt,$L/ovl.gp2_sch.txt,pd_gp2,rls" "$L/lava_pd_rls_gp2.csv"

# ---- A1-AFR (g1000_afr) ----
mk_info "$L/info.afr.txt" pd_afr 1488 196430 "$M/pd_rizig_afr.sumstats.gz" rls_afr 2176 153313 "$M/rls_afr.sumstats.gz"; mk_ovl "$L/ovl.afr.txt" pd_afr rls_afr
runscan AFR g1000_afr "afr,$L/info.afr.txt,$L/ovl.afr.txt,pd_afr,rls_afr" "$L/lava_pd_rls_afr.csv"

# ---- A1-AMR (g1000_amr) ----
mk_info "$L/info.amr.txt" pd_amr 740 741 "$M/pd_loesch_amr.sumstats.gz" rls_amr 2024 91902 "$M/rls_amr.sumstats.gz"; mk_ovl "$L/ovl.amr.txt" pd_amr rls_amr
runscan AMR g1000_amr "amr,$L/info.amr.txt,$L/ovl.amr.txt,pd_amr,rls_amr" "$L/lava_pd_rls_amr.csv"

echo "### ALL DONE — results: results/lava/lava_pd_rls_{gp2,afr,amr}.csv ; tell Claude to finalize the cross-ancestry story ###"
