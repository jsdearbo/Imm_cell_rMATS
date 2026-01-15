#!/bin/bash
set -euo pipefail

PLOTTER="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/plot_rmats_volcano.py"
RESULT_ROOT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results"
DEFAULT_COMPARISONS=(
  "B_vs_T"
  "B_vs_Macrophage"
  "T_vs_Macrophage"
)
HIGHLIGHTS=()
EXCLUDES=()
TOP_N=10
OUTPUT_SUBDIR="plots"
PREFIX=""
DPI=160
FORCE=0
comparisons=()

while (( $# > 0 )); do
  case "$1" in
    -c|--comparison)
      [[ $# -ge 2 ]] || { echo "Error: --comparison needs an argument" >&2; exit 1; }
      comparisons+=("$2")
      shift 2
      ;;
    -H|--highlight)
      [[ $# -ge 2 ]] || { echo "Error: --highlight needs an argument" >&2; exit 1; }
      HIGHLIGHTS+=("$2")
      shift 2
      ;;
    -x|--exclude)
      [[ $# -ge 2 ]] || { echo "Error: --exclude needs an argument" >&2; exit 1; }
      EXCLUDES+=("$2")
      shift 2
      ;;
    --top-n)
      TOP_N="$2"
      shift 2
      ;;
    --output-subdir)
      OUTPUT_SUBDIR="$2"
      shift 2
      ;;
    --prefix)
      PREFIX="$2"
      shift 2
      ;;
    --dpi)
      DPI="$2"
      shift 2
      ;;
    -f|--force)
      FORCE=1
      shift
      ;;
    -h|--help)
      echo "Usage: $0 [options]"
      echo "  -c, --comparison <label>   Add comparison (repeatable)."
      echo "  -H, --highlight <gene>     Highlight gene (repeatable)."
      echo "  -x, --exclude <gene>       Exclude gene (repeatable)."
      echo "      --top-n <int>           Top hits to label (default: 10)."
      echo "      --output-subdir <name>  Subdir for plots (default: plots)."
      echo "      --prefix <text>         Filename prefix."
      echo "      --dpi <int>             Figure DPI (default: 160)."
      echo "  -f, --force                Regenerate plots even if sentinel exists."
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

[[ -f "${PLOTTER}" ]] || { echo "Error: plot script missing at ${PLOTTER}" >&2; exit 1; }

for comp in "${comparisons[@]}"; do
  comp_dir="${RESULT_ROOT}/${comp}"
  if [[ ! -d "${comp_dir}" ]]; then
    echo "[warn] missing comparison directory: ${comp_dir}" >&2
    continue
  fi
  sentinel="${comp_dir}/${OUTPUT_SUBDIR}/plots_complete.flag"
  if [[ -f "${sentinel}" && ${FORCE} -eq 0 ]]; then
    echo "=== Skipping ${comp} plots (sentinel present) ==="
    continue
  fi
  mkdir -p "${comp_dir}/${OUTPUT_SUBDIR}"
  cmd=(
    python "${PLOTTER}" "${comp_dir}"
    --top-n "${TOP_N}"
    --output-subdir "${OUTPUT_SUBDIR}"
    --dpi "${DPI}"
  )
  if [[ -n "${PREFIX}" ]]; then
    cmd+=(--prefix "${PREFIX}")
  fi
  for gene in "${HIGHLIGHTS[@]}"; do
    cmd+=(--highlight "${gene}")
  done
  for gene in "${EXCLUDES[@]}"; do
    cmd+=(--exclude "${gene}")
  done
  echo "=== Plotting ${comp} ==="
  "${cmd[@]}"
  touch "${sentinel}"
  echo "=== Completed ${comp} ==="
  echo
  sleep 1
done
