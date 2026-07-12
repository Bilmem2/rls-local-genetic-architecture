#!/usr/bin/env bash
HERE="$(cd "$(dirname "$0")" && pwd)"  # this script's directory
cd /path/to/analysis || exit 1
G="${1:-rlssig}"
export PAIR_GROUP="$G"
Rscript "$HERE/panel_lava2.R" > "results/_logs/panel_lava_${G}.log" 2>&1
echo "PANEL_LAVA_EXIT $? [$G]" >> "results/_logs/panel_lava_${G}.log"
tail -6 "results/_logs/panel_lava_${G}.log"
