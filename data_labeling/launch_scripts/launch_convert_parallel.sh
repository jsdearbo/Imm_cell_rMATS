#!/usr/bin/env bash
set -euo pipefail

# --------------------
# Configuration
# --------------------
PYTHON_SCRIPT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/data_conversion_fast.py"
GENOME_FASTA="/mnt/data_sda/jdearborn/genomes/mouse/mm10/fasta/GRCm38.primary_assembly.genome.fa"
MANIFEST="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/outputs_soft_per_sample/soft/manifest.tsv"
TRANSCRIPTS_ALL="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/stringtie_out/dominant_isoforms/all_celltypes_transcripts.tsv"
BASE_OUT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/results/labelled_data/windowed_multitask_fast"

mkdir -p "${BASE_OUT}"
cp "${MANIFEST}" "${BASE_OUT}/manifest.tsv"

# Function to process a single chromosome (or set of chromosomes)
run_chrom() {
    local SPLIT=$1
    local CHROM=$2
    
    # Write to a temp subdirectory to avoid shard name collisions (shard_0000.npz)
    local SUB_DIR="${BASE_OUT}/${SPLIT}/${CHROM}"
    
    echo "[${SPLIT}] Starting ${CHROM}..."
    
    python "${PYTHON_SCRIPT}" \
      --fasta "${GENOME_FASTA}" \
      --transcripts "${TRANSCRIPTS_ALL}" \
      --manifest "${MANIFEST}" \
      --out-dir "${BASE_OUT}" \
      --split "${SPLIT}/${CHROM}" \
      --chroms "${CHROM}" \
      --window-bp 524288 \
      --subset-mode balanced \
      --shard-size 256 > "${BASE_OUT}/log_${CHROM}.txt" 2>&1
      
    echo "[${SPLIT}] Finished ${CHROM}"
}

export -f run_chrom
export PYTHON_SCRIPT GENOME_FASTA MANIFEST TRANSCRIPTS_ALL BASE_OUT

# --------------------
# EXECUTION
# --------------------

# List of chromosomes
CHROMS_TRAIN="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY"

echo "Launching Parallel Jobs for TRAIN..."

# Method 1: If you have GNU parallel installed (Recommended)
# parallel -j 5 run_chrom train ::: ${CHROMS_TRAIN}

# Method 2: Simple bash background loops (If no GNU parallel)
MAX_JOBS=5
count=0
for chrom in ${CHROMS_TRAIN}; do
    run_chrom train "$chrom" &
    
    # Simple job throttling
    count=$((count+1))
    if [ "$count" -ge "$MAX_JOBS" ]; then
        wait -n # Wait for at least one job to finish
        count=$((count-1))
    fi
done
wait # Wait for all remaining

# --------------------
# MERGING
# --------------------
# Move all shards from train/chrX/*.npz to train/*.npz
echo "Merging shards..."
mkdir -p "${BASE_OUT}/train"
find "${BASE_OUT}/train" -mindepth 2 -name "*.npz" -print0 | xargs -0 -I {} mv --backup=t {} "${BASE_OUT}/train/"

# Clean up empty chrom dirs
find "${BASE_OUT}/train" -type d -empty -delete

echo "Done."