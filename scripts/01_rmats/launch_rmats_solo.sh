#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

RMATS_FILE="/mnt/netfiles/installations/rmats_turbo_v4_3_0/rmats.py"
GTF_FILE="/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf"
SEQ_TYPE="paired"
READ_LENGTH=151
NTHREAD=8

OUTPUT_ROOT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results"
TMP_ROOT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/rmats/tmp_single"
FORCE=0
groups=()

declare -A GROUP_FILES=(
  [B]="${REPO_ROOT}/config/B_Cell_bams.txt"
  [T]="${REPO_ROOT}/config/T_Cell_bams.txt"
  [Macrophage]="${REPO_ROOT}/config/Macrophage_bams.txt"
)

DEFAULT_GROUPS=(B T Macrophage)

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
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

run_single_rmats() {
  local group="$1"
  local manifest="${GROUP_FILES[${group}]}"
  [[ -f "${manifest}" ]] || die "Manifest missing: ${manifest}"

  local out_dir="${OUTPUT_ROOT}/${group}_only"
  local run_tmp="${TMP_ROOT}/${group}_$(date +%Y%m%d_%H%M%S)"
  local sentinel="${out_dir}/rmats_complete.flag"

  if [[ -f "${sentinel}" && ${FORCE} -eq 0 ]]; then
    echo "=== Skipping ${group} (sentinel present) ==="
    return
  fi

  if [[ ${FORCE} -eq 1 ]]; then
    echo "=== Cleaning old results for ${group} ==="
    rm -rf "${out_dir}"
    rm -f "${sentinel}"
  fi

  mkdir -p "${out_dir}" "${run_tmp}"

  echo "=== Running ${group} only ==="
  python "${RMATS_FILE}" \
    --b1 "${manifest}" \
    --gtf "${GTF_FILE}" \
    -t "${SEQ_TYPE}" \
    --allow-clipping \
    --variable-read-length \
    --readLength "${READ_LENGTH}" \
    --nthread "${NTHREAD}" \
    --od "${out_dir}" \
    --tmp "${run_tmp}"

  touch "${sentinel}"
  echo "=== Completed ${group} only ==="
}

while (( $# > 0 )); do
  case "$1" in
    -f|--force)
      FORCE=1
      ;;
    *)
      groups+=("$(parse_group "$1")")
      ;;
  esac
  shift
done

if (( ${#groups[@]} == 0 )); then
  groups=("${DEFAULT_GROUPS[@]}")
fi

for grp in "${groups[@]}"; do
  run_single_rmats "${grp}"
  echo
  sleep 2
done