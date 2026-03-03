#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# StringTie per-replicate pipeline
#
# Workflow (best-practice for scaffold GTF construction):
#   1. Run StringTie -e per individual BAM  (reference-only quantification)
#   2. Merge per-cell-type GTFs with stringtie --merge
#   3. Re-quantify each replicate against the merged model (-e -B)
#
# The final merged GTFs are what downstream dominant_isoform_labler.py consumes.
###############################################################################

########################
# USER SETTINGS
########################
GUIDE_GTF="/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf"
BAM_LIST_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt"
OUTDIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_per_rep"
THREADS=12
PARALLEL_JOBS=${PARALLEL_JOBS:-12}

# Library is unstranded
STRAND_OPT=""

CELL_TYPES=(B_Cell T_Cell Macrophage)
declare -A BAM_LISTS=(
  [B_Cell]="${BAM_LIST_DIR}/B_Cell_bams.txt"
  [T_Cell]="${BAM_LIST_DIR}/T_Cell_bams.txt"
  [Macrophage]="${BAM_LIST_DIR}/Macrophage_bams.txt"
)

########################
# Derived paths
########################
PER_REP_DIR="${OUTDIR}/per_replicate"
MERGED_DIR="${OUTDIR}/merged"
REQUANT_DIR="${OUTDIR}/requantified"

mkdir -p "${PER_REP_DIR}" "${MERGED_DIR}" "${REQUANT_DIR}"

if [[ "${PARALLEL_JOBS}" -lt 1 ]]; then PARALLEL_JOBS=1; fi

########################
# Helper: parse comma-delimited BAM list file into an array
########################
parse_bam_list() {
  local list_file=$1
  tr ',' '\n' < "${list_file}" \
    | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' \
    | sed '/^$/d'
}

########################
# Helper: run a function with GNU parallel-style batching
########################
wait_batch() {
  local -n _pids=$1
  for pid in "${_pids[@]}"; do
    wait "${pid}" || echo "[WARN] PID ${pid} exited non-zero" >&2
  done
  _pids=()
}

########################
# STEP 1: StringTie per replicate (-e, reference-only)
########################
echo "================================================================"
echo "[STEP 1] Running StringTie per replicate (reference-only, -e)"
echo "================================================================"

job_pids=()
job_count=0

for cell in "${CELL_TYPES[@]}"; do
  list_file="${BAM_LISTS[${cell}]}"
  if [[ ! -f "${list_file}" ]]; then
    echo "[ERROR] BAM list not found: ${list_file}" >&2
    continue
  fi

  while IFS= read -r bam; do
    sample=$(basename "${bam}" .bam)
    out_gtf="${PER_REP_DIR}/${cell}_${sample}.gtf"

    if [[ ! -f "${bam}" ]]; then
      echo "[WARN] BAM missing, skipping: ${bam}" >&2
      continue
    fi

    # Index if needed
    if [[ ! -f "${bam}.bai" ]]; then
      echo "[INFO] Indexing ${bam}..."
      samtools index -@ 2 "${bam}"
    fi

    (
      echo "[INFO] StringTie -e: ${cell} / ${sample}"
      stringtie "${bam}" \
        -G "${GUIDE_GTF}" \
        -e \
        -o "${out_gtf}" \
        -p 2 \
        ${STRAND_OPT} \
        -l "${cell}_${sample}"
    ) &
    job_pids+=("$!")
    job_count=$((job_count + 1))

    if [[ "${job_count}" -ge "${PARALLEL_JOBS}" ]]; then
      wait_batch job_pids
      job_count=0
    fi
  done < <(parse_bam_list "${list_file}")
done

wait_batch job_pids
echo "[STEP 1] Done."

########################
# STEP 2: StringTie --merge per cell type
########################
echo "================================================================"
echo "[STEP 2] Merging per-replicate GTFs (stringtie --merge)"
echo "================================================================"

for cell in "${CELL_TYPES[@]}"; do
  gtf_list="${MERGED_DIR}/${cell}_gtf_list.txt"
  merged_gtf="${MERGED_DIR}/${cell}_merged.gtf"

  # Collect all per-replicate GTFs for this cell type
  find "${PER_REP_DIR}" -name "${cell}_*.gtf" -type f | sort > "${gtf_list}"

  n_gtfs=$(wc -l < "${gtf_list}")
  if [[ "${n_gtfs}" -eq 0 ]]; then
    echo "[WARN] No per-replicate GTFs found for ${cell}, skipping merge." >&2
    continue
  fi

  echo "[INFO] Merging ${n_gtfs} GTFs for ${cell}..."
  stringtie --merge \
    -G "${GUIDE_GTF}" \
    -o "${merged_gtf}" \
    -p "${THREADS}" \
    ${STRAND_OPT} \
    "${gtf_list}"

  echo "[INFO] Wrote merged model: ${merged_gtf}"
done

echo "[STEP 2] Done."

########################
# STEP 3: Re-quantify each replicate against the merged model
########################
echo "================================================================"
echo "[STEP 3] Re-quantifying replicates against merged models (-e)"
echo "================================================================"

job_pids=()
job_count=0

for cell in "${CELL_TYPES[@]}"; do
  merged_gtf="${MERGED_DIR}/${cell}_merged.gtf"
  if [[ ! -f "${merged_gtf}" ]]; then
    echo "[WARN] No merged GTF for ${cell}, skipping re-quantification." >&2
    continue
  fi

  list_file="${BAM_LISTS[${cell}]}"

  while IFS= read -r bam; do
    sample=$(basename "${bam}" .bam)
    rq_dir="${REQUANT_DIR}/${cell}_${sample}"
    mkdir -p "${rq_dir}"

    (
      echo "[INFO] Re-quantify: ${cell} / ${sample}"
      stringtie "${bam}" \
        -G "${merged_gtf}" \
        -e \
        -B \
        -o "${rq_dir}/${cell}_${sample}.gtf" \
        -A "${rq_dir}/${cell}_${sample}_gene_abund.tsv" \
        -p 2 \
        ${STRAND_OPT} \
        -l "${cell}_${sample}"
    ) &
    job_pids+=("$!")
    job_count=$((job_count + 1))

    if [[ "${job_count}" -ge "${PARALLEL_JOBS}" ]]; then
      wait_batch job_pids
      job_count=0
    fi
  done < <(parse_bam_list "${list_file}")
done

wait_batch job_pids
echo "[STEP 3] Done."

########################
# STEP 4: Build one aggregate GTF per cell type (average TPM across replicates)
########################
echo "================================================================"
echo "[STEP 4] Aggregating re-quantified GTFs per cell type"
echo "================================================================"

# For each cell type, pick the merged model GTF (from step 2) as the scaffold.
# The per-replicate re-quantified GTFs (step 3) carry sample-level TPM/cov.
# We create a simple symlink so downstream dominant_isoform_labler.py can
# consume ${OUTDIR}/${cell}_aggregate.gtf in the same way it did before.
#
# If you need mean-TPM across replicates, use the Python script (dominant_isoform_labler.py)
# pointed at the requantified directory, or use prepDE.py from StringTie.

for cell in "${CELL_TYPES[@]}"; do
  merged_gtf="${MERGED_DIR}/${cell}_merged.gtf"
  agg_link="${OUTDIR}/${cell}_aggregate.gtf"

  if [[ -f "${merged_gtf}" ]]; then
    ln -sf "${merged_gtf}" "${agg_link}"
    echo "[INFO] Linked ${agg_link} -> ${merged_gtf}"
  fi
done

echo "[STEP 4] Done."
echo "================================================================"
echo "[DONE] Per-replicate StringTie pipeline complete."
echo "  Per-replicate GTFs : ${PER_REP_DIR}/"
echo "  Merged models      : ${MERGED_DIR}/"
echo "  Re-quantified      : ${REQUANT_DIR}/"
echo "  Aggregate symlinks : ${OUTDIR}/*_aggregate.gtf"
echo "================================================================"
