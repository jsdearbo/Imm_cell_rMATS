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
# 1) consume merged BAM per cell type and run StringTie
########################
CELL_TYPES=(B_Cell T_Cell Macrophage)

mkdir -p "${OUTDIR}"

run_stringtie() {
  local prefix=$1
  local merged_bam="${BAM_DIR}/${prefix}_merged.bam"
  local out_gtf="${OUTDIR}/${prefix}_aggregate.gtf"

  if [[ ! -f "${merged_bam}" ]]; then
    echo "[ERROR] Expected merged BAM missing: ${merged_bam}" >&2
    return 1
  fi

  if [[ ! -f "${merged_bam}.bai" ]]; then
    echo "[INFO] Indexing ${merged_bam}..."
    samtools index -@ "${THREADS}" "${merged_bam}"
  fi

  echo "[INFO] Running StringTie on ${prefix} merged BAM..."
  stringtie "${merged_bam}" \
    -G "${GUIDE_GTF}" \
    -o "${out_gtf}" \
    -p "${THREADS}" \
    ${STRAND_OPT} \
    -l "${prefix}_agg"
}

for cell in "${CELL_TYPES[@]}"; do
  run_stringtie "${cell}"
done

echo "[INFO] Done."