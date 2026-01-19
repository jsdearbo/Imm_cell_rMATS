#!/usr/bin/env bash
set -euo pipefail

PYTHON_SCRIPT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/data_labeling/data_conversion.py"

GENOME_FASTA="/mnt/netfiles/genomes/STAR/mm10/GRCm38.primary_assembly.genome.fa"

MANIFEST_SRC="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/data_labeling/outputs_soft_per_sample/soft/manifest.tsv"

BASE_OUT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/labelled_data/windowed"

# Chromosome splits
TRAIN_CHROMS=$(echo chr{1..9},chr{12..19},chrX,chrY | tr ' ' ',')
VAL_CHROMS="chr10"
TEST_CHROMS="chr11"

declare -A TRANSCRIPTS=(
  [B_Cell]="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/B_Cell_transcripts.tsv"
  [T_Cell]="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/T_Cell_transcripts.tsv"
  [Macrophage]="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/Macrophage_transcripts.tsv"
)

for CELL in "${!TRANSCRIPTS[@]}"; do
  echo "=========================================="
  echo "[INFO] Converting dataset for ${CELL}"
  echo "=========================================="

  OUT_DIR="${BASE_OUT}/${CELL}"
  mkdir -p "${OUT_DIR}"

  # Copy manifest once per dataset
  cp "${MANIFEST_SRC}" "${OUT_DIR}/manifest.tsv"

  # --------------------
  # TRAIN
  # --------------------
  python "${PYTHON_SCRIPT}" \
    --fasta "${GENOME_FASTA}" \
    --transcripts "${TRANSCRIPTS[$CELL]}" \
    --out-dir "${OUT_DIR}" \
    --split train \
    --chroms "${TRAIN_CHROMS}" \
    --window-bp 524288 \
    --bin-size 32 \
    --stride-bp 196608 \
    --shard-size 256

  # --------------------
  # VAL
  # --------------------
  python "${PYTHON_SCRIPT}" \
    --fasta "${GENOME_FASTA}" \
    --transcripts "${TRANSCRIPTS[$CELL]}" \
    --out-dir "${OUT_DIR}" \
    --split val \
    --chroms "${VAL_CHROMS}" \
    --window-bp 524288 \
    --bin-size 32 \
    --stride-bp 196608 \
    --shard-size 128

  # --------------------
  # TEST
  # --------------------
  python "${PYTHON_SCRIPT}" \
    --fasta "${GENOME_FASTA}" \
    --transcripts "${TRANSCRIPTS[$CELL]}" \
    --out-dir "${OUT_DIR}" \
    --split test \
    --chroms "${TEST_CHROMS}" \
    --window-bp 524288 \
    --bin-size 32 \
    --stride-bp 196608 \
    --shard-size 128

done

echo "[INFO] Full data conversion complete (train/val/test)."
