#!/usr/bin/env bash
# make_figure1.sh — generate all Figure 1 panels
# Outputs go to analysis/figure1/figure1_out/
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/figure1_out"
mkdir -p "${OUT_DIR}"

run_panel() {
    local name="$1"
    local script="$2"
    echo ""
    echo "──────────────────────────────────────────"
    echo "  Panel ${name}"
    echo "──────────────────────────────────────────"
    python3 "${script}"
}

run_panel "A — RNA-seq cohort overview"        "${SCRIPT_DIR}/panel_A_cohort.py"
run_panel "B — StringTie binary skeleton"      "${SCRIPT_DIR}/panel_B_skeleton.py"
run_panel "C — rMATS PSI integration"          "${SCRIPT_DIR}/panel_C_integration.py"
run_panel "D — Example locus"                  "${SCRIPT_DIR}/panel_D_locus_example.py"
run_panel "E — PSI label distribution"         "${SCRIPT_DIR}/panel_E_distribution.py"
run_panel "F — Event counts by class"          "${SCRIPT_DIR}/panel_F_event_counts.py"

echo ""
echo "══════════════════════════════════════════"
echo "  All panels complete. PDFs in:"
echo "  ${OUT_DIR}"
ls "${OUT_DIR}"/*.pdf 2>/dev/null || echo "  (no PDFs found — check errors above)"
