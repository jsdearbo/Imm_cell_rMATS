#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SOFT_SCRIPT="${PROJECT_ROOT}/src/add_confidence_to_bed.py"

# Input Roots
RMATS_ROOT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/results"
STRINGTIE_ROOT="${RMATS_ROOT}/stringtie_out/dominant_isoforms"

# Output Root (separate from old tests)
OUTPUT_ROOT="${SCRIPT_DIR}/outputs_soft_per_sample_confidence"
SOFT_OUT_DIR="${OUTPUT_ROOT}/soft"

mkdir -p "${SOFT_OUT_DIR}"

declare -A CELL_TO_SUBDIR=(
  [B_Cell]="B_only"
  [T_Cell]="T_only"
  [Macrophage]="Macrophage_only"
)

declare -A GTF_PATHS=(
  [B_Cell]="${STRINGTIE_ROOT}/B_Cell_dominant.gtf"
  [T_Cell]="${STRINGTIE_ROOT}/T_Cell_dominant.gtf"
  [Macrophage]="${STRINGTIE_ROOT}/Macrophage_dominant.gtf"
)

declare -A EVENT_FILES=(
  [SE]="SE.MATS.JC.txt"
  [RI]="RI.MATS.JC.txt"
  [A3SS]="A3SS.MATS.JC.txt"
  [A5SS]="A5SS.MATS.JC.txt"
)

for CT in B_Cell T_Cell Macrophage; do
  echo "=========================================="
  echo "[INFO] Testing per-sample soft labels for ${CT}"
  echo "=========================================="

  subdir="${CELL_TO_SUBDIR[${CT}]:-}"
  if [[ -z "${subdir}" ]]; then
    echo "[WARN] No rMATS directory mapping for ${CT}; skipping"
    continue
  fi

  rmats_dir="${RMATS_ROOT}/${subdir}"
  sentinel="${rmats_dir}/rmats_complete.flag"

  if [[ ! -f "${sentinel}" ]]; then
    echo "[WARN] rMATS run incomplete for ${CT} (missing ${sentinel}); skipping"
    continue
  fi

  gtf_path="${GTF_PATHS[${CT}]}"
  if [[ ! -f "${gtf_path}" ]]; then
    echo "[WARN] Missing GTF for ${CT}: ${gtf_path}; skipping"
    continue
  fi

  se_path="${rmats_dir}/${EVENT_FILES[SE]}"
  ri_path="${rmats_dir}/${EVENT_FILES[RI]}"
  a3_path="${rmats_dir}/${EVENT_FILES[A3SS]}"
  a5_path="${rmats_dir}/${EVENT_FILES[A5SS]}"

  if [[ ! -f "${se_path}" ]]; then
    echo "[WARN] Missing SE data for ${CT}: ${se_path}; skipping"
    continue
  fi

  for event in RI A3SS A5SS; do
    path="${rmats_dir}/${EVENT_FILES[${event}]}"
    if [[ ! -f "${path}" ]]; then
      echo "[WARN] Missing ${event} data for ${CT}: ${path}"
    fi
  done

  # Output prefix (directory + cell type stem)
  out_prefix="${SOFT_OUT_DIR}/${CT}_soft"

  echo "[INFO] Emitting per-sample soft-label BEDs for ${CT}"
  python "${SOFT_SCRIPT}" \
    --rmats-se "${se_path}" \
    --rmats-ri "${ri_path}" \
    --rmats-a3ss "${a3_path}" \
    --rmats-a5ss "${a5_path}" \
    --gtf "${gtf_path}" \
    --out-prefix "${out_prefix}"

done

echo "[INFO] Per-sample soft label test complete."
