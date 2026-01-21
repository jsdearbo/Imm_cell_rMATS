#!/usr/bin/env bash
set -euo pipefail

PYTHON_SCRIPT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/make_manifest.py"
SOFT_LABEL_ROOT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/outputs_soft_per_sample/soft"

python "$PYTHON_SCRIPT" \
  --soft-label-root "$SOFT_LABEL_ROOT" \
  --out "$SOFT_LABEL_ROOT/manifest.tsv"

echo "Manifest completed successfully."