#!/usr/bin/env bash
set -euo pipefail

########################
# USER SETTINGS
########################

GUIDE_GTF="/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf"
BAM_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/merge_bams"
OUTDIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out"
THREADS=12

# Anti-artifact controls (can be overridden via env vars)
STRICT_REFERENCE=${STRICT_REFERENCE:-1}   # 1: only quantify reference transcripts (-e)
MIN_COVERAGE=${MIN_COVERAGE:-2.0}         # -c
MIN_JUNCTION_COV=${MIN_JUNCTION_COV:-2}   # -j
MIN_ISOFORM_FRAC=${MIN_ISOFORM_FRAC:-0.05} # -f
MIN_TRANSCRIPT_LEN=${MIN_TRANSCRIPT_LEN:-200} # -m

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
  local out_abund="${OUTDIR}/${prefix}_gene_abundance.tsv"

  if [[ ! -f "${merged_bam}" ]]; then
    echo "[ERROR] Expected merged BAM missing: ${merged_bam}" >&2
    return 1
  fi

  if [[ ! -f "${merged_bam}.bai" ]]; then
    echo "[INFO] Indexing ${merged_bam}..."
    samtools index -@ "${THREADS}" "${merged_bam}"
  fi

  echo "[INFO] Running StringTie on ${prefix} merged BAM..."
  if [[ "${STRICT_REFERENCE}" == "1" ]]; then
    stringtie "${merged_bam}" \
      -G "${GUIDE_GTF}" \
      -e \
      -o "${out_gtf}" \
      -A "${out_abund}" \
      -p "${THREADS}" \
      ${STRAND_OPT} \
      -l "${prefix}_agg"
  else
    stringtie "${merged_bam}" \
      -G "${GUIDE_GTF}" \
      -o "${out_gtf}" \
      -A "${out_abund}" \
      -p "${THREADS}" \
      -c "${MIN_COVERAGE}" \
      -j "${MIN_JUNCTION_COV}" \
      -f "${MIN_ISOFORM_FRAC}" \
      -m "${MIN_TRANSCRIPT_LEN}" \
      ${STRAND_OPT} \
      -l "${prefix}_agg"
  fi
}

for cell in "${CELL_TYPES[@]}"; do
  run_stringtie "${cell}"
done

echo "[INFO] Done."