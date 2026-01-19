#!/usr/bin/env bash
set -euo pipefail

PYTHON_SCRIPT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/data_labeling/data_conversion.py"

GENOME_FASTA="/mnt/netfiles/genomes/STAR/mm10/GRCm38.primary_assembly.genome.fa"

MANIFEST_SRC="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/data_labeling/outputs_soft_per_sample/soft/manifest.tsv"

BASE_OUT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/labelled_data/windowed"

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

  # Reuse the same manifest, filtered implicitly by transcripts + task index
  cp "${MANIFEST_SRC}" "${OUT_DIR}/manifest.tsv"

  python "${PYTHON_SCRIPT}" \
    --fasta "${GENOME_FASTA}" \
    --transcripts "${TRANSCRIPTS[$CELL]}" \
    --out-dir "${OUT_DIR}" \
    --split train \
    --window-bp 196608 \
    --bin-size 32 \
    --stride-bp 98304 \
    --shard-size 256

done

echo "[INFO] Full data conversion complete."
