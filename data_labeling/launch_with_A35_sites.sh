#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Scripts
HARD_SCRIPT="${SCRIPT_DIR}/test_hard_labeling.py"
SOFT_SCRIPT="${SCRIPT_DIR}/test_soft_labeling.py"
VALID_SCRIPT="${SCRIPT_DIR}/validation_tool.py"

# Input Roots
RMATS_ROOT="${PROJECT_ROOT}/results"
STRINGTIE_ROOT="${RMATS_ROOT}/stringtie_out/dominant_isoforms"
OUTPUT_ROOT="${SCRIPT_DIR}/outputs"
HARD_OUT_DIR="${OUTPUT_ROOT}/hard"
SOFT_OUT_DIR="${OUTPUT_ROOT}/soft"
VALID_OUT_DIR="${OUTPUT_ROOT}/validation"

mkdir -p "${HARD_OUT_DIR}" "${SOFT_OUT_DIR}" "${VALID_OUT_DIR}"

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
  echo "[INFO] Processing ${CT}"
  echo "=========================================="

  subdir="${CELL_TO_SUBDIR[${CT}]}"
  rmats_dir="${RMATS_ROOT}/${subdir}"
  sentinel="${rmats_dir}/rmats_complete.flag"

  if [[ -z "${subdir}" ]]; then
    echo "[WARN] No rMATS directory mapping for ${CT}; skipping"
    continue
  fi

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

  for event in SE RI A3SS A5SS; do
    path="${rmats_dir}/${EVENT_FILES[${event}]}"
    if [[ ! -f "${path}" ]]; then
      echo "[WARN] Missing ${event} data for ${CT}: ${path}"
    fi
  done

  hard_bed="${HARD_OUT_DIR}/${CT}_hard_labels.bed"
  hard_summary="${HARD_OUT_DIR}/${CT}_hard_label_summary.tsv"
  soft_bed="${SOFT_OUT_DIR}/${CT}_soft_labels.bed"
  soft_summary="${SOFT_OUT_DIR}/${CT}_soft_label_summary.tsv"
  validation_bed="${VALID_OUT_DIR}/${CT}_validation_diff.bed"

  # 1. Generate Hard Labels
  python "${HARD_SCRIPT}" \
    --rmats-se "${se_path}" \
    --rmats-ri "${ri_path}" \
    --rmats-a3ss "${a3_path}" \
    --rmats-a5ss "${a5_path}" \
    --gtf "${gtf_path}" \
    --out-bed "${hard_bed}" \
    --summary "${hard_summary}"

  # 2. Generate Soft Labels
  python "${SOFT_SCRIPT}" \
    --rmats-se "${se_path}" \
    --rmats-ri "${ri_path}" \
    --rmats-a3ss "${a3_path}" \
    --rmats-a5ss "${a5_path}" \
    --gtf "${gtf_path}" \
    --out-bed "${soft_bed}" \
    --summary "${soft_summary}"

  # 3. Generate Validation Disagreements
  echo "[INFO] Generating Validation Track..."
  python "${VALID_SCRIPT}" \
    --gtf "${gtf_path}" \
    --hard-bed "${hard_bed}" \
    --out-bed "${validation_bed}"

done

echo "[INFO] All cell types complete."