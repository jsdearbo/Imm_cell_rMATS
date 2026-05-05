#!/usr/bin/env bash
# ==============================================================================
# pipeline_config.sh - Central configuration for the splicing label pipeline.
#
# Source this from any launch script or from run_full_pipeline.sh.
# All paths and tunable parameters live here so they aren't scattered
# across multiple shell scripts.
# ==============================================================================

# ---- Root directories ----
PROJECT_ROOT="/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS"
RESULTS_ROOT="${PROJECT_ROOT}/results"
LABELING_ROOT="${PROJECT_ROOT}/data_labeling"

# ---- Genome ----
GENOME_FASTA="/mnt/data_sda/jdearborn/genomes/mouse/mm10/fasta/GRCm38.primary_assembly.genome.fa"
GUIDE_GTF="/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf"

# ---- Cell types (comma-separated for Python, space-separated array for bash) ----
CELL_TYPES_CSV="B_Cell,T_Cell,Macrophage"
CELL_TYPES_ARR=(B_Cell T_Cell Macrophage)

# ---- StringTie / dominant isoform outputs ----
STRINGTIE_DIR="${RESULTS_ROOT}/stringtie_per_rep"
DOMINANT_GTF_DIR="${STRINGTIE_DIR}/dominant_isoforms"

# ---- Step 1: GTF -> TSV conversion ----
GTF_TO_TSV_SCRIPT="${LABELING_ROOT}/gtf_to_tsv.py"
TRANSCRIPTS_ALL="${DOMINANT_GTF_DIR}/all_celltypes_transcripts.tsv"

# ---- Step 2: Soft labeling (rMATS PSI -> BED) ----
SOFT_LABEL_SCRIPT="${LABELING_ROOT}/test_soft_labeling.py"
RMATS_ROOT="${RESULTS_ROOT}"
GTF_ROOT="${DOMINANT_GTF_DIR}"
SOFT_LABEL_OUT_DIR="${LABELING_ROOT}/results/PSI_labelled_bed"
PER_SAMPLE_OUT_DIR="${SOFT_LABEL_OUT_DIR}/per_sample"
MANIFEST_OUT="${LABELING_ROOT}/launch_scripts/manifest.generated.tsv"

# ---- Step 3: NPZ shard conversion ----
CONVERT_SCRIPT="${LABELING_ROOT}/data_conversion_fast.py"
NPZ_OUT_DIR="${LABELING_ROOT}/results/labelled_data/multitask_balanced"
WINDOW_BP=524288
SHARD_SIZE=256
JITTER_BP=16
SUBSET_MODE="balanced"

# ---- Chromosome splits ----
CHROMS_TRAIN="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY"
CHROMS_VAL="chr10"
CHROMS_TEST="chr11"

# ---- Parallelism ----
MAX_PARALLEL_CHROMS=5
SHARD_PREFIX_STEP=10000
