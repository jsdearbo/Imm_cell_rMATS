#!/usr/bin/env bash
set -euo pipefail

# --------------------
# Configuration
# --------------------
PYTHON_SCRIPT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/data_conversion_fast.py"
GENOME_FASTA="/mnt/data_sda/jdearborn/genomes/mouse/mm10/fasta/GRCm38.primary_assembly.genome.fa"
MANIFEST="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/outputs_soft_per_sample/soft/manifest.tsv"
TRANSCRIPTS_ALL="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/stringtie_out/dominant_isoforms/all_celltypes_transcripts.tsv"

# Separate folder for debug
BASE_OUT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/labelled_data/DEBUG_multitask_robust"

mkdir -p "${BASE_OUT}"
cp "${MANIFEST}" "${BASE_OUT}/manifest.tsv"

echo "Starting DEBUG Generation (Small subset)..."

# --------------------
# DEBUG RUN
# --------------------
# 1. Runs only on chr19 (usually smaller/denser gene content)
# 2. Uses --max-examples 50 to stop quickly
# 3. Uses 'balanced' mode to ensure we see both signal and empty windows

python "${PYTHON_SCRIPT}" \
  --fasta "${GENOME_FASTA}" \
  --transcripts "${TRANSCRIPTS_ALL}" \
  --manifest "${MANIFEST}" \
  --out-dir "${BASE_OUT}" \
  --split train \
  --chroms "chr19" \
  --window-bp 524288 \
  --subset-mode balanced \
  --shard-size 10 \
  --jitter-bp 16 \
  --augment-k 1 \
  --max-examples 50

# Also generate a tiny VAL set for the trainer to use
python "${PYTHON_SCRIPT}" \
  --fasta "${GENOME_FASTA}" \
  --transcripts "${TRANSCRIPTS_ALL}" \
  --manifest "${MANIFEST}" \
  --out-dir "${BASE_OUT}" \
  --split val \
  --chroms "chr19" \
  --window-bp 524288 \
  --subset-mode balanced \
  --shard-size 10 \
  --jitter-bp 16 \
  --max-examples 20

echo "[INFO] Debug generation complete. Check folder: ${BASE_OUT}"