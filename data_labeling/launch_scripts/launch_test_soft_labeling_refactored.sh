#!/usr/bin/env bash
set -euo pipefail

# inputs
PIPELINE_SCRIPT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/test_soft_labeling.py"
RMATS_ROOT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/results"
GTF_ROOT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/results/stringtie_per_rep/dominant_isoforms"

# outputs
OUT_DIR="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/PSI_labelled_bed"
PER_SAMPLE_OUT_DIR="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/PSI_labelled_bed/per_sample"
MANIFEST_OUT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/launch_scripts/manifest.generated.tsv"

CELL_TYPES="B_Cell,T_Cell,Macrophage"

mkdir -p "${OUT_DIR}" "${PER_SAMPLE_OUT_DIR}"

echo "=========================================="
echo "[INFO] Running refactored soft-label pipeline"
echo "[INFO] 1) Regenerate per-sample BEDs from rMATS"
echo "[INFO] 2) Generate manifest"
echo "[INFO] 3) Build master BEDs"
echo "=========================================="

python "${PIPELINE_SCRIPT}" \
  --out-dir "${OUT_DIR}" \
  --regenerate-beds \
  --rmats-root "${RMATS_ROOT}" \
  --gtf-root "${GTF_ROOT}" \
  --per-sample-out-dir "${PER_SAMPLE_OUT_DIR}" \
  --manifest-out "${MANIFEST_OUT}" \
  --cell-types "${CELL_TYPES}"

echo "[INFO] Pipeline complete."
echo "[INFO] Manifest: ${MANIFEST_OUT}"
echo "[INFO] Master BEDs: ${OUT_DIR}/master_beds"
