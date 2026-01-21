#!/usr/bin/env bash
set -euo pipefail

# --------------------
# Configuration
# --------------------
PYTHON_SCRIPT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/data_conversion.py"
GENOME_FASTA="/mnt/data_sda/jdearborn/genomes/mouse/mm10/fasta/GRCm38.primary_assembly.genome.fa"
MANIFEST="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/outputs_soft_per_sample/soft/manifest.tsv"

# IMPORTANT: transcript union TSV
# You should create this once by concatenating and de-duplicating
TRANSCRIPTS_ALL="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/stringtie_out/dominant_isoforms/all_celltypes_transcripts.tsv"

BASE_OUT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/labelled_data/windowed_multitask"

# --------------------
# Chromosome splits
# --------------------
TRAIN_CHROMS=$(echo chr{1..9},chr{12..19},chrX,chrY | tr ' ' ',')
VAL_CHROMS="chr10"
TEST_CHROMS="chr11"

mkdir -p "${BASE_OUT}"
cp "${MANIFEST}" "${BASE_OUT}/manifest.tsv"

# --------------------
# TRAIN
# --------------------
python "${PYTHON_SCRIPT}" \
  --fasta "${GENOME_FASTA}" \
  --transcripts "${TRANSCRIPTS_ALL}" \
  --out-dir "${BASE_OUT}" \
  --split train \
  --chroms "${TRAIN_CHROMS}" \
  --window-bp 524288 \
  --bin-size 32 \
  --stride-bp 262144 \
  --shard-size 256 \
  --manifest "${MANIFEST}"

# --------------------
# VAL
# --------------------
python "${PYTHON_SCRIPT}" \
  --fasta "${GENOME_FASTA}" \
  --transcripts "${TRANSCRIPTS_ALL}" \
  --out-dir "${BASE_OUT}" \
  --split val \
  --chroms "${VAL_CHROMS}" \
  --window-bp 524288 \
  --bin-size 32 \
  --stride-bp 262144 \
  --shard-size 128 \
  --manifest "${MANIFEST}"

# --------------------
# TEST
# --------------------
python "${PYTHON_SCRIPT}" \
  --fasta "${GENOME_FASTA}" \
  --transcripts "${TRANSCRIPTS_ALL}" \
  --out-dir "${BASE_OUT}" \
  --split test \
  --chroms "${TEST_CHROMS}" \
  --window-bp 524288 \
  --bin-size 32 \
  --stride-bp 262144 \
  --shard-size 128 \
  --manifest "${MANIFEST}"

echo "[INFO] Multitask data conversion complete (train/val/test)."
