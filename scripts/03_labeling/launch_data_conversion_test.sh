#!/usr/bin/env bash
set -euo pipefail

PYTHON_SCRIPT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/data_labeling/data_conversion.py"

# Genome FASTA (indexed with .fai)
GENOME_FASTA="/mnt/netfiles/genomes/STAR/mm10/GRCm38.primary_assembly.genome.fa"

# Transcripts TSV (ONE cell type for smoke test)
TRANSCRIPTS_TSV="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/B_Cell_transcripts.tsv"

# Output directory
OUT_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/labelled_data/windowed_smoketest/B_Cell"

# Manifest listing per-sample BEDs
MANIFEST_SRC="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/data_labeling/outputs_soft_per_sample/soft/manifest.tsv"

mkdir -p "${OUT_DIR}"
cp "${MANIFEST_SRC}" "${OUT_DIR}/manifest.tsv"

python "${PYTHON_SCRIPT}" \
  --fasta "${GENOME_FASTA}" \
  --transcripts "${TRANSCRIPTS_TSV}" \
  --out-dir "${OUT_DIR}" \
  --split train \
  --window-bp 131072 \
  --bin-size 32 \
  --stride-bp 65536 \
  --shard-size 16 \
  --max-examples 32

echo "[INFO] Smoke test conversion complete."
