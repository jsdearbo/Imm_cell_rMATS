#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# run_full_pipeline.sh - Unified splicing-label pipeline
#
# Runs the data-labeling steps in sequence:
#
#   Step 1:  Convert dominant-isoform GTFs to transcript TSVs
#   Step 2:  Generate per-sample soft-label BED files from rMATS PSI values,
#            then aggregate into master BEDs and write a manifest
#   Step 3:  Convert BED labels + genome FASTA into NPZ shards for training
#
# Upstream prerequisites (NOT run by this script):
#   - StringTie per-replicate annotation   (make_gtf/stringtie_per_replicate.sh)
#   - Dominant isoform selection            (make_gtf/dominant_isoform_labler.py)
#   - rMATS differential splicing analysis
#
# Usage:
#   bash run_full_pipeline.sh                # run all three steps
#   bash run_full_pipeline.sh --from 2       # skip GTF-to-TSV, start at labeling
#   bash run_full_pipeline.sh --only 3       # run only NPZ conversion
#   bash run_full_pipeline.sh --dry-run      # print what would run, then exit
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG="${REPO_ROOT}/config/pipeline_config.sh"

if [[ ! -f "${CONFIG}" ]]; then
    echo "[ERROR] Missing config: ${CONFIG}" >&2
    exit 1
fi
# shellcheck source=pipeline_config.sh
source "${CONFIG}"

# ---------- argument parsing ----------
FROM_STEP=1
TO_STEP=3
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)    FROM_STEP="$2"; shift 2 ;;
        --to)      TO_STEP="$2";   shift 2 ;;
        --only)    FROM_STEP="$2"; TO_STEP="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help)
            sed -n '3,/^# ====/{ /^# ====/d; s/^# //; s/^#//; p; }' "$0"
            exit 0 ;;
        *) echo "[ERROR] Unknown option: $1" >&2; exit 1 ;;
    esac
done

should_run() { [[ "$1" -ge "${FROM_STEP}" && "$1" -le "${TO_STEP}" ]]; }

banner() {
    echo ""
    echo "============================================================"
    echo "  [STEP $1]  $2"
    echo "============================================================"
}

# ---------- pre-flight checks ----------
preflight_ok=true

check_file() {
    if [[ ! -f "$1" ]]; then
        echo "[PREFLIGHT] Missing: $1" >&2
        preflight_ok=false
    fi
}

check_file "${GTF_TO_TSV_SCRIPT}"
check_file "${SOFT_LABEL_SCRIPT}"
check_file "${CONVERT_SCRIPT}"

if should_run 1; then
    for ct in "${CELL_TYPES_ARR[@]}"; do
        check_file "${DOMINANT_GTF_DIR}/${ct}_dominant.gtf"
    done
fi

if should_run 3; then
    check_file "${GENOME_FASTA}"
fi

if [[ "${preflight_ok}" != "true" ]]; then
    echo "[ERROR] Pre-flight checks failed. Resolve missing files above." >&2
    exit 1
fi

if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "[DRY RUN] Would execute steps ${FROM_STEP}–${TO_STEP}. Exiting."
    exit 0
fi

PIPELINE_START=$(date +%s)

# ======================================================================
# STEP 1: Convert dominant-isoform GTFs to transcript TSVs
# ======================================================================
if should_run 1; then
    banner 1 "GTF -> TSV conversion"

    # Remove stale union file so it builds fresh
    rm -f "${TRANSCRIPTS_ALL}"

    for ct in "${CELL_TYPES_ARR[@]}"; do
        GTF_IN="${DOMINANT_GTF_DIR}/${ct}_dominant.gtf"
        TSV_OUT="${DOMINANT_GTF_DIR}/${ct}_transcripts.tsv"

        echo "[Step 1] Converting ${ct}: ${GTF_IN}"
        python "${GTF_TO_TSV_SCRIPT}" \
            --gtf "${GTF_IN}" \
            --out "${TSV_OUT}" \
            --union-out "${TRANSCRIPTS_ALL}"
    done

    echo "[Step 1] Done.  Union TSV: ${TRANSCRIPTS_ALL}"
fi

# ======================================================================
# STEP 2: Generate soft-label BED files from rMATS PSI values
# ======================================================================
if should_run 2; then
    banner 2 "Soft labeling (rMATS PSI -> per-sample BEDs + master BEDs)"

    mkdir -p "${SOFT_LABEL_OUT_DIR}" "${PER_SAMPLE_OUT_DIR}"

    python "${SOFT_LABEL_SCRIPT}" \
        --out-dir "${SOFT_LABEL_OUT_DIR}" \
        --regenerate-beds \
        --rmats-root "${RMATS_ROOT}" \
        --gtf-root "${GTF_ROOT}" \
        --per-sample-out-dir "${PER_SAMPLE_OUT_DIR}" \
        --manifest-out "${MANIFEST_OUT}" \
        --cell-types "${CELL_TYPES_CSV}"

    echo "[Step 2] Done.  Manifest: ${MANIFEST_OUT}"
fi

# ======================================================================
# STEP 3: Convert BED labels + FASTA -> NPZ shards  (parallel by chrom)
# ======================================================================
if should_run 3; then
    banner 3 "NPZ shard conversion (parallel by chromosome)"

    # Use the manifest produced by step 2 (or a pre-existing one)
    MANIFEST="${MANIFEST_OUT}"
    if [[ ! -f "${MANIFEST}" ]]; then
        echo "[ERROR] Manifest not found: ${MANIFEST}" >&2
        echo "        Run step 2 first, or point MANIFEST_OUT at an existing file." >&2
        exit 1
    fi

    mkdir -p "${NPZ_OUT_DIR}"
    cp "${MANIFEST}" "${NPZ_OUT_DIR}/manifest.tsv"

    # --- helper launched in subshells ---
    run_chrom() {
        local SPLIT=$1  CHROM=$2  PREFIX=$3
        echo "[${SPLIT}] Starting ${CHROM} (shard prefix ${PREFIX})..."

        python "${CONVERT_SCRIPT}" \
            --fasta     "${GENOME_FASTA}" \
            --transcripts "${TRANSCRIPTS_ALL}" \
            --manifest  "${MANIFEST}" \
            --out-dir   "${NPZ_OUT_DIR}" \
            --split     "${SPLIT}/${CHROM}" \
            --chroms    "${CHROM}" \
            --window-bp "${WINDOW_BP}" \
            --subset-mode "${SUBSET_MODE}" \
            --shard-size  "${SHARD_SIZE}" \
            --jitter-bp   "${JITTER_BP}" \
            --shard-prefix "${PREFIX}" \
            > "${NPZ_OUT_DIR}/log_${CHROM}.txt" 2>&1

        echo "[${SPLIT}] Finished ${CHROM}"
    }
    export -f run_chrom
    export CONVERT_SCRIPT GENOME_FASTA TRANSCRIPTS_ALL NPZ_OUT_DIR
    export WINDOW_BP SUBSET_MODE SHARD_SIZE JITTER_BP

    # --- TRAIN chromosomes (throttled) ---
    echo "[Step 3] Launching TRAIN jobs..."
    PREFIX=0
    count=0
    for chrom in ${CHROMS_TRAIN}; do
        run_chrom train "${chrom}" "${PREFIX}" &
        PREFIX=$((PREFIX + SHARD_PREFIX_STEP))
        count=$((count + 1))
        if [[ "${count}" -ge "${MAX_PARALLEL_CHROMS}" ]]; then
            wait -n
            count=$((count - 1))
        fi
    done

    # --- VAL + TEST ---
    run_chrom val  "${CHROMS_VAL}"  0 &
    run_chrom test "${CHROMS_TEST}" 0 &

    wait  # all background jobs

    # --- Merge per-chrom subdirs into flat split dirs ---
    echo "[Step 3] Merging shards..."
    for SPLIT in train val test; do
        mkdir -p "${NPZ_OUT_DIR}/${SPLIT}"
        find "${NPZ_OUT_DIR}/${SPLIT}" -mindepth 2 -name "shard_*.npz" -print0 \
            | xargs -0 -r -I {} mv {} "${NPZ_OUT_DIR}/${SPLIT}/"
        find "${NPZ_OUT_DIR}/${SPLIT}" -type d -empty -delete 2>/dev/null || true
    done

    echo "[Step 3] Done.  Shards in ${NPZ_OUT_DIR}/{train,val,test}/"
fi

# ---------- summary ----------
PIPELINE_END=$(date +%s)
ELAPSED=$((PIPELINE_END - PIPELINE_START))
MINS=$((ELAPSED / 60))
SECS=$((ELAPSED % 60))

echo ""
echo "============================================================"
echo "  Pipeline complete  (steps ${FROM_STEP}–${TO_STEP})"
echo "  Wall time: ${MINS}m ${SECS}s"
echo "============================================================"
