#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS"
PYTHON_SCRIPT="${BASE_DIR}/src/gtf_to_tsv.py"

B_CELL_GTF="${BASE_DIR}/results/stringtie_per_rep/dominant_isoforms/B_Cell_dominant.gtf"
B_CELL_TSV="${BASE_DIR}/results/stringtie_per_rep/dominant_isoforms/B_Cell_transcripts.tsv"

T_CELL_GTF="${BASE_DIR}/results/stringtie_per_rep/dominant_isoforms/T_Cell_dominant.gtf"
T_CELL_TSV="${BASE_DIR}/results/stringtie_per_rep/dominant_isoforms/T_Cell_transcripts.tsv"

MAC_GTF="${BASE_DIR}/results/stringtie_per_rep/dominant_isoforms/Macrophage_dominant.gtf"
MAC_TSV="${BASE_DIR}/results/stringtie_per_rep/dominant_isoforms/Macrophage_transcripts.tsv"

UNIFIED_TSV="${BASE_DIR}/results/stringtie_per_rep/dominant_isoforms/all_celltypes_transcripts.tsv"

convert_jobs=(
  "B Cell|$B_CELL_GTF|$B_CELL_TSV"
  "T Cell|$T_CELL_GTF|$T_CELL_TSV"
  "Macrophage|$MAC_GTF|$MAC_TSV"
)

echo "Starting GTF to TSV conversions (total ${#convert_jobs[@]} jobs)..."

for job in "${convert_jobs[@]}"; do
  IFS='|' read -r label gtf_path tsv_path <<< "$job"
  echo "[$label] Converting $gtf_path"
  python "$PYTHON_SCRIPT" --gtf "$gtf_path" --out "$tsv_path" --union-out "$UNIFIED_TSV"
  echo "[$label] Wrote $tsv_path"
done

echo "All conversions completed successfully."