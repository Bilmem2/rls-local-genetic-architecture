#!/usr/bin/env bash
HERE="$(cd "$(dirname "$0")" && pwd)"  # this script's directory
cd /path/to/analysis || exit 1
mkdir -p results/_logs
Rscript "$HERE/panel_lava.R" > results/_logs/panel_lava.log 2>&1
echo "PANEL_LAVA_EXIT $?" >> results/_logs/panel_lava.log
echo "=== tail ==="; tail -30 results/_logs/panel_lava.log
