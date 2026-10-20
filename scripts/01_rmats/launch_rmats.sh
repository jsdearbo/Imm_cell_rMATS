#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Core rMATS configuration shared by every comparison run.
RMATS_FILE="/mnt/netfiles/installations/rmats_turbo_v4_3_0/rmats.py"
GTF_FILE="/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf"
SEQ_TYPE="paired"

# Updated to 151 based on your "Mode" read length check.
READ_LENGTH=151

NTHREAD=8
OUTPUT_ROOT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results"
TMP_ROOT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/rmats/tmp"
FORCE=0

comparisons=()

while (( $# > 0 )); do
  case "$1" in
    -f|--force)
      FORCE=1
      shift
      ;;
    *)
      comparisons+=("$1")
      shift
      ;;
  esac
done

# Central registry of sample manifests (config/ in repo root)
declare -A GROUP_FILES=(
  [B]="${REPO_ROOT}/config/B_Cell_bams.txt"
  [T]="${REPO_ROOT}/config/T_Cell_bams.txt"
  [Macrophage]="${REPO_ROOT}/config/Macrophage_bams.txt"
)

DEFAULT_COMPARISONS=(
  "B_vs_T"
  "B_vs_Macrophage"
  "T_vs_Macrophage"
)

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
  export LD_LIBRARY_PATH
fi

mkdir -p "${OUTPUT_ROOT}" "${TMP_ROOT}"

die() {
  echo "Error: $*" >&2
  exit 1
}

parse_group() {
  local label="$1"
  local upper="${label^}"
  case "${upper}" in
    B) echo "B" ;;
    T) echo "T" ;;
    Macrophage|Mac) echo "Macrophage" ;;
    *) die "Unknown group '${label}'. Available: ${!GROUP_FILES[*]}" ;;
  esac
}

run_rmats() {
  local group1="$1"
  local group2="$2"
  local comp_label="${group1}_vs_${group2}"
  local b1="${GROUP_FILES[${group1}]}"
  local b2="${GROUP_FILES[${group2}]}"
  [[ -f "${b1}" ]] || die "Manifest missing: ${b1}"
  [[ -f "${b2}" ]] || die "Manifest missing: ${b2}"
  local out_dir="${OUTPUT_ROOT}/${comp_label}"
  local run_tmp="${TMP_ROOT}/${comp_label}_$(date +%Y%m%d_%H%M%S)"
  local sentinel="${out_dir}/rmats_complete.flag"

  if [[ -f "${sentinel}" && ${FORCE} -eq 0 ]]; then
    echo "=== Skipping ${comp_label} (sentinel present) ==="
    return
  fi
  
  # If forcing, clean up the old directory to avoid mixing results
  if [[ ${FORCE} -eq 1 ]]; then
    echo "=== Cleaning old results for ${comp_label} ==="
    rm -rf "${out_dir}"
    rm -f "${sentinel}"
  fi
  
  mkdir -p "${out_dir}" "${run_tmp}"

  echo "=== Running ${comp_label} ==="
  python "${RMATS_FILE}" \
    --b1 "${b1}" \
    --b2 "${b2}" \
    --gtf "${GTF_FILE}" \
    -t "${SEQ_TYPE}" \
    --allow-clipping \
    --variable-read-length \
    --readLength "${READ_LENGTH}" \
    --nthread "${NTHREAD}" \
    --od "${out_dir}" \
    --tmp "${run_tmp}"
  
  touch "${sentinel}"
  echo "=== Completed ${comp_label} ==="
}

if (( ${#comparisons[@]} == 0 )); then
  comparisons=("${DEFAULT_COMPARISONS[@]}")
fi

for comparison in "${comparisons[@]}"; do
  if [[ "${comparison}" != *_vs_* ]]; then
    die "Comparison '${comparison}' must follow <group1>_vs_<group2> syntax"
  fi
  group1_raw="${comparison%%_vs_*}"
  group2_raw="${comparison##*_vs_}"
  group1=$(parse_group "${group1_raw}")
  group2=$(parse_group "${group2_raw}")
  run_rmats "${group1}" "${group2}"
  echo
  sleep 2
done