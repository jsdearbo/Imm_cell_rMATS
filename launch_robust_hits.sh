#!/bin/bash
set -euo pipefail

EXTRACTOR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/extract_robust_hits.py"
RESULT_ROOT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results"
DEFAULT_COMPARISONS=(
  "B_vs_T"
  "B_vs_Macrophage"
  "T_vs_Macrophage"
)
FDR=0.0000000005
DPSI=0.1
TOP=50
OUTPUT_NAME="robust_events_summary.csv"
FORCE=0
comparisons=()

while (( $# > 0 )); do
  case "$1" in
    -c|--comparison)
      [[ $# -ge 2 ]] || { echo "Error: --comparison needs an argument" >&2; exit 1; }
      comparisons+=("$2")
      shift 2
      ;;
    --fdr)
      FDR="$2"
      shift 2
      ;;
    --dpsi)
      DPSI="$2"
      shift 2
      ;;
    --top)
      TOP="$2"
      shift 2
      ;;
    -o|--output)
      OUTPUT_NAME="$2"
      shift 2
      ;;
    -f|--force)
      FORCE=1
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [options]"
      echo "  -c, --comparison <label>   Add comparison (repeatable)."
      echo "      --fdr <float>          FDR threshold (default: 0.05)."
      echo "      --dpsi <float>         |dPSI| threshold (default: 0.1)."
      echo "      --top <int>            Rows to show in summary (default: 50)."
      echo "  -o, --output <filename>   Override output CSV name."
      echo "  -f, --force               Regenerate even if sentinel exists."
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if (( ${#comparisons[@]} == 0 )); then
  comparisons=("${DEFAULT_COMPARISONS[@]}")
fi

[[ -f "${EXTRACTOR}" ]] || { echo "Error: extractor script missing at ${EXTRACTOR}" >&2; exit 1; }

for comp in "${comparisons[@]}"; do
  comp_dir="${RESULT_ROOT}/${comp}"
  if [[ ! -d "${comp_dir}" ]]; then
    echo "[warn] missing comparison directory: ${comp_dir}" >&2
    continue
  fi
  sentinel="${comp_dir}/robust_complete.flag"
  if [[ -f "${sentinel}" && ${FORCE} -eq 0 ]]; then
    echo "=== Skipping ${comp} robust extraction (sentinel present) ==="
    continue
  fi
  cmd=(
    python "${EXTRACTOR}" "${comp_dir}"
    --fdr "${FDR}"
    --dpsi "${DPSI}"
    --top "${TOP}"
    --output "${OUTPUT_NAME}"
  )
  echo "=== Extracting robust events for ${comp} ==="
  "${cmd[@]}"
  touch "${sentinel}"
  echo "=== Completed ${comp} ==="
  echo
  sleep 1
done
