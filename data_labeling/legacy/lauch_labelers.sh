#!/usr/bin/env bash
set -euo pipefail

BASE="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq"
TEST_DIR="${BASE}/compute_psi/hard_and_soft_label_testing"

HARD_SCRIPT="${TEST_DIR}/test_hard_labeling.py"
SOFT_SCRIPT="${TEST_DIR}/test_soft_labeling.py"

# Define Paths for SE (Skipped Exon)
declare -A SE_PATHS=(
  [B_Cell]="${BASE}/rmats/B_only/SE.MATS.JC.txt"
  [T_Cell]="${BASE}/rmats/T_only/SE.MATS.JC.txt"
  [Macrophage]="${BASE}/rmats/Macrophage_only/SE.MATS.JC.txt"
)

# Define Paths for RI (Retained Intron)
declare -A RI_PATHS=(
  [B_Cell]="${BASE}/rmats/B_only/RI.MATS.JC.txt"
  [T_Cell]="${BASE}/rmats/T_only/RI.MATS.JC.txt"
  [Macrophage]="${BASE}/rmats/Macrophage_only/RI.MATS.JC.txt"
)

declare -A GTF_PATHS=(
  [B_Cell]="${BASE}/compute_psi/stringtie_out/dominant_isoforms/B_Cell_dominant.gtf"
  [T_Cell]="${BASE}/compute_psi/stringtie_out/dominant_isoforms/T_Cell_dominant.gtf"
  [Macrophage]="${BASE}/compute_psi/stringtie_out/dominant_isoforms/Macrophage_dominant.gtf"
)

for CT in B_Cell T_Cell Macrophage; do
  echo "[INFO] Processing ${CT}"

  HARD_OUT_DIR="${TEST_DIR}/hard_output_test"
  SOFT_OUT_DIR="${TEST_DIR}/soft_output_test"
  mkdir -p "${HARD_OUT_DIR}" "${SOFT_OUT_DIR}"

  # HARD LABELS
  python "${HARD_SCRIPT}" \
    --rmats-se "${SE_PATHS[$CT]}" \
    --rmats-ri "${RI_PATHS[$CT]}" \
    --gtf "${GTF_PATHS[$CT]}" \
    --out-bed "${HARD_OUT_DIR}/${CT}_hard_labels.bed" \
    --summary "${HARD_OUT_DIR}/${CT}_hard_label_summary.tsv"

  # SOFT LABELS
  python "${SOFT_SCRIPT}" \
    --rmats-se "${SE_PATHS[$CT]}" \
    --rmats-ri "${RI_PATHS[$CT]}" \
    --gtf "${GTF_PATHS[$CT]}" \
    --out-bed "${SOFT_OUT_DIR}/${CT}_soft_labels.bed" \
    --summary "${SOFT_OUT_DIR}/${CT}_soft_label_summary.tsv"
done

echo "[INFO] All cell types complete."