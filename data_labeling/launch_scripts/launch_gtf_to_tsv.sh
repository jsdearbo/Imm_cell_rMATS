#!/usr/bin/env bash
set -euo pipefail

PYTHON_SCRIPT="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/data_labeling/gtf_to_tsv.py"

B_CELL_GTF="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/B_Cell_dominant.gtf"
B_CELL_TSV="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/B_Cell_transcripts.tsv"

T_CELL_GTF="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/T_Cell_dominant.gtf"
T_CELL_TSV="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/T_Cell_transcripts.tsv"

MAC_GTF="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/Macrophage_dominant.gtf"
MAC_TSV="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/Macrophage_transcripts.tsv"

UNIFIED_TSV="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out/dominant_isoforms/all_celltypes_transcripts.tsv"

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