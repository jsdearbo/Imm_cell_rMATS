#!/usr/bin/env bash
set -euo pipefail

# --------------------
# Configuration
# --------------------
# Update this to point to the new FIXED python script
PYTHON_SCRIPT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/src/data_conversion.py"

GENOME_FASTA="/mnt/data_sda/jdearborn/genomes/mouse/mm10/fasta/GRCm38.primary_assembly.genome.fa"
MANIFEST="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/outputs_soft_per_sample/soft/manifest.tsv"
TRANSCRIPTS_ALL="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/stringtie_out/dominant_isoforms/all_celltypes_transcripts.tsv"

# Output directory for the REAL run
BASE_OUT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/labelled_data/windowed_multitask_balanced"

# --------------------
# Chromosome splits
# --------------------
TRAIN_CHROMS=$(echo chr{1..9},chr{12..19},chrX,chrY | tr ' ' ',')
VAL_CHROMS="chr10"
TEST_CHROMS="chr11"

mkdir -p "${BASE_OUT}"
cp "${MANIFEST}" "${BASE_OUT}/manifest.tsv"

echo "Starting Full Dataset Generation..."
echo "Output: ${BASE_OUT}"

# --------------------
# TRAIN
# --------------------
echo "--- Processing TRAIN ---"
python "${PYTHON_SCRIPT}" \
  --fasta "${GENOME_FASTA}" \
  --transcripts "${TRANSCRIPTS_ALL}" \
  --manifest "${MANIFEST}" \
  --out-dir "${BASE_OUT}" \
  --split train \
  --chroms "${TRAIN_CHROMS}" \
  --window-bp 524288 \
  --subset-mode balanced \
  --shard-size 256

# --------------------
# VAL
# --------------------
echo "--- Processing VAL ---"
python "${PYTHON_SCRIPT}" \
  --fasta "${GENOME_FASTA}" \
  --transcripts "${TRANSCRIPTS_ALL}" \
  --manifest "${MANIFEST}" \
  --out-dir "${BASE_OUT}" \
  --split val \
  --chroms "${VAL_CHROMS}" \
  --window-bp 524288 \
  --subset-mode balanced \
  --shard-size 128

# --------------------
# TEST
# --------------------
echo "--- Processing TEST ---"
python "${PYTHON_SCRIPT}" \
  --fasta "${GENOME_FASTA}" \
  --transcripts "${TRANSCRIPTS_ALL}" \
  --manifest "${MANIFEST}" \
  --out-dir "${BASE_OUT}" \
  --split test \
  --chroms "${TEST_CHROMS}" \
  --window-bp 524288 \
  --subset-mode balanced \
  --shard-size 128

echo "[INFO] Multitask data conversion complete."