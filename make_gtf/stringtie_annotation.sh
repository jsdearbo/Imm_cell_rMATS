#!/usr/bin/env bash
set -euo pipefail

########################
# USER SETTINGS
########################

GUIDE_GTF="/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf"
BAM_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/merge_bams"
OUTDIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out"
THREADS=16

# strandedness
#STRAND_OPT="--rf"   # fr-firststrand
#STRAND_OPT="--fr"   # fr-secondstrand
STRAND_OPT=""        # unstranded

########################
# 1) merged BAM per cell type, then StringTie on merged BAM
########################
CELL_TYPES=(B_Cell T_Cell Macrophage)

mkdir -p "${OUTDIR}"
shopt -s nullglob

merge_and_stringtie() {
  local prefix="$1"
  local merged_bam="${OUTDIR}/${prefix}_merged.bam"
  local out_gtf="${OUTDIR}/${prefix}_aggregate.gtf"

  echo "[INFO] Merging ${prefix} BAMs..."
  samtools merge -@ "${THREADS}" "${merged_bam}" "${BAM_DIR}/${prefix}"*.bam
  samtools index "${merged_bam}"

  echo "[INFO] Running StringTie on ${prefix} merged BAM..."
  stringtie "${merged_bam}" \
    -G "${GUIDE_GTF}" \
    -o "${out_gtf}" \
    -p "${THREADS}" \
    ${STRAND_OPT} \
    -l "${prefix}_agg"
}

for cell in "${CELL_TYPES[@]}"; do
  merge_and_stringtie "${cell}"
done

echo "[INFO] Done."