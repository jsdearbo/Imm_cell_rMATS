#!/usr/bin/env bash
set -euo pipefail

########################
# USER SETTINGS
########################

GUIDE_GTF="/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf"
BAM_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/merge_bams"
OUTDIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out"
THREADS=12
PARALLEL_JOBS=${PARALLEL_JOBS:-3}

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

if [[ "${PARALLEL_JOBS}" -lt 1 ]]; then
  PARALLEL_JOBS=1
fi
if [[ "${PARALLEL_JOBS}" -gt "${#CELL_TYPES[@]}" ]]; then
  PARALLEL_JOBS=${#CELL_TYPES[@]}
fi

THREADS_PER_JOB=$(( THREADS / PARALLEL_JOBS ))
if [[ "${THREADS_PER_JOB}" -lt 1 ]]; then
  THREADS_PER_JOB=1
fi

echo "[INFO] THREADS=${THREADS}, PARALLEL_JOBS=${PARALLEL_JOBS}, THREADS_PER_JOB=${THREADS_PER_JOB}"

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
      -p "${THREADS_PER_JOB}" \
      ${STRAND_OPT} \
      -l "${prefix}_agg"
  else
    stringtie "${merged_bam}" \
      -G "${GUIDE_GTF}" \
      -o "${out_gtf}" \
      -A "${out_abund}" \
      -p "${THREADS_PER_JOB}" \
      -c "${MIN_COVERAGE}" \
      -j "${MIN_JUNCTION_COV}" \
      -f "${MIN_ISOFORM_FRAC}" \
      -m "${MIN_TRANSCRIPT_LEN}" \
      ${STRAND_OPT} \
      -l "${prefix}_agg"
  fi
}

job_count=0
pids=()
for cell in "${CELL_TYPES[@]}"; do
  run_stringtie "${cell}" &
  pids+=("$!")
  job_count=$((job_count + 1))

  if [[ "${job_count}" -ge "${PARALLEL_JOBS}" ]]; then
    for pid in "${pids[@]}"; do
      wait "${pid}"
    done
    pids=()
    job_count=0
  fi
done

for pid in "${pids[@]}"; do
  wait "${pid}"
done

echo "[INFO] Done."