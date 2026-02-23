#!/usr/bin/env bash
set -euo pipefail

# --------------------
# Configuration
# --------------------
PYTHON_SCRIPT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/data_conversion_fast.py"
GENOME_FASTA="/mnt/data_sda/jdearborn/genomes/mouse/mm10/fasta/GRCm38.primary_assembly.genome.fa"
MANIFEST="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/outputs_soft_per_sample/soft/manifest.tsv"
TRANSCRIPTS_ALL="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/stringtie_out/dominant_isoforms/all_celltypes_transcripts.tsv"
BASE_OUT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/labelled_data/multitask_confidence_event_only_augment_4"

mkdir -p "${BASE_OUT}"
cp "${MANIFEST}" "${BASE_OUT}/manifest.tsv"

# Function: process one chromosome
run_chrom() {
    local SPLIT=$1
    local CHROM=$2
    local PREFIX=$3
    
    local SUB_DIR="${BASE_OUT}/${SPLIT}/${CHROM}"
    
    echo "[${SPLIT}] Starting ${CHROM} (Shard Prefix: ${PREFIX})..."
    
    # NOTE: Set --augment-k > 1 here to generate multiple shifted copies per window
    # OR reduce --stride-bp to increase overlap.
    
    python "${PYTHON_SCRIPT}" \
      --fasta "${GENOME_FASTA}" \
      --transcripts "${TRANSCRIPTS_ALL}" \
      --manifest "${MANIFEST}" \
      --out-dir "${BASE_OUT}" \
      --split "${SPLIT}/${CHROM}" \
      --chroms "${CHROM}" \
      --window-bp 524288 \
      --subset-mode event_only \
      --shard-size 256 \
      --jitter-bp 16 \
      --augment-k 4 \
      --shard-prefix "${PREFIX}" > "${BASE_OUT}/log_${CHROM}.txt" 2>&1
      
    echo "[${SPLIT}] Finished ${CHROM}"
}

export -f run_chrom
export PYTHON_SCRIPT GENOME_FASTA MANIFEST TRANSCRIPTS_ALL BASE_OUT

# --------------------
# EXECUTION
# --------------------

# 1. TRAIN SPLIT (Parallel)
CHROMS_TRAIN="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY"

echo "Launching Parallel Jobs for TRAIN..."

# Launch jobs in background
# We increment prefix by 10000 for each chromosome (Increased to handle potential augmentation data scale)
PREFIX=0
STEP=10000
MAX_JOBS=5
count=0

for chrom in ${CHROMS_TRAIN}; do
    run_chrom train "$chrom" "$PREFIX" &
    
    PREFIX=$((PREFIX + STEP))
    count=$((count+1))
    
    # Simple throttling
    if [ "$count" -ge "$MAX_JOBS" ]; then
        wait -n
        count=$((count-1))
    fi
done

# 2. VAL SPLIT (Single Job - chr10)
# Run in background but part of the pool
run_chrom val "chr10" 0 &
count=$((count+1))

# 3. TEST SPLIT (Single Job - chr11)
run_chrom test "chr11" 0 &

wait

# --------------------
# MERGING (Clean)
# --------------------
# Now we can just move the files because names are unique
echo "Merging shards..."
mkdir -p "${BASE_OUT}/train" "${BASE_OUT}/val" "${BASE_OUT}/test"

# Function to merge a split
merge_split() {
    local SPLIT=$1
    echo "  Merging ${SPLIT}..."
    find "${BASE_OUT}/${SPLIT}" -mindepth 2 -name "shard_*.npz" -print0 | xargs -0 -r -I {} mv {} "${BASE_OUT}/${SPLIT}/"
    find "${BASE_OUT}/${SPLIT}" -type d -empty -delete
}

merge_split "train"
merge_split "val"
merge_split "test"

echo "Done."