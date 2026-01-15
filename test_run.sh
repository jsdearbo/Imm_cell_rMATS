#!/bin/bash

# Define variables
RMATS_FILE="/mnt/netfiles/installations/rmats_turbo_v4_3_0/rmats.py"
B1_FILE="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/B_Cell_bams.txt"
# B1_FILE="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/T_Cell_bams.txt"
# B1_FILE="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/Macrophage_bams.txt"

# B2_FILE="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/B_Cell_bams.txt"  # use B2 if running comparisons
# B2_FILE="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/T_Cell_bams.txt"
B2_FILE="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/Macrophage_bams.txt"

GTF_FILE="/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf"
READ_LENGTH=150
NTHREAD=8
# OUTPUT_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/B_only"
# OUTPUT_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/T_only"
# OUTPUT_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/Macrophage_only"
# OUTPUT_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/B_vs_T"
OUTPUT_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/B_vs_Macrophage"
# OUTPUT_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/T_vs_Macrophage"

TMP_DIR="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/rmats/tmp"
SEQ_TYPE="paired"

# Set the library path
# export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

# Run rmats.py with the specified parameters-- add allow-clipping and variable-read-length options
python "$RMATS_FILE" \
    --b1 "$B1_FILE" \
    --b2 "$B2_FILE" \
    --gtf "$GTF_FILE" \
    -t "$SEQ_TYPE" \
    --allow-clipping \
    --variable-read-length \
    --readLength "$READ_LENGTH" \
    --nthread "$NTHREAD" \
    --od "$OUTPUT_DIR" \
    --tmp "$TMP_DIR"
