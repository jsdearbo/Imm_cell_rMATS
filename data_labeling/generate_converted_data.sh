#!/bin/bash
set -euo pipefail

BASE="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq"
TEST_DIR="${BASE}/compute_psi/hard_and_soft_label_testing"
REF_GTF="/bobross/STAR_Metadata/mm10/gencode.vM23.annotation.gtf"

# Define Input BEDs (From previous validation step)
B_HARD="${TEST_DIR}/hard_output_test/B_Cell_hard_labels.bed"
T_HARD="${TEST_DIR}/hard_output_test/T_Cell_hard_labels.bed"
M_HARD="${TEST_DIR}/hard_output_test/Macrophage_hard_labels.bed"

B_SOFT="${TEST_DIR}/soft_output_test/B_Cell_soft_labels.bed"
T_SOFT="${TEST_DIR}/soft_output_test/T_Cell_soft_labels.bed"
M_SOFT="${TEST_DIR}/soft_output_test/Macrophage_soft_labels.bed"

# 1. Generate HARD Target TSV
echo "Generating Hard Label TSV..."
python /bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/compute_psi/hard_and_soft_label_testing/dataset_conversion.py\
  --ref-gtf "$REF_GTF" \
  --mode hard \
  --b-bed "$B_HARD" \
  --t-bed "$T_HARD" \
  --m-bed "$M_HARD" \
  --out-tsv "${TEST_DIR}/training_targets_hard.shift0.tsv"

# 2. Generate SOFT Target TSV
echo "Generating Soft Label TSV..."
python /bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/compute_psi/hard_and_soft_label_testing/dataset_conversion.py \
  --ref-gtf "$REF_GTF" \
  --mode soft \
  --b-bed "$B_SOFT" \
  --t-bed "$T_SOFT" \
  --m-bed "$M_SOFT" \
  --out-tsv "${TEST_DIR}/training_targets_soft.shift0.tsv"

echo "Done. You can now point MultiTaskClassificationDataset to:"
echo "  Hard: ${TEST_DIR}/training_targets_hard"
echo "  Soft: ${TEST_DIR}/training_targets_soft"