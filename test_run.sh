#!/bin/bash

# Define variables
RMATS_FILE="/bobross/jdearborn/install_tools/rmats_turbo_v4_3_0/rmats.py"
B1_FILE="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/txt/B_Cell_tracks.txt"
#B1_FILE="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/txt/T_Cell_tracks.txt"
#B1_FILE="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/txt/Macrophage_tracks.txt"

#B2_FILE="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/txt/B_Cell_tracks.txt"  # use B2 if running comparisons
#B2_FILE="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/txt/T_Cell_tracks.txt"
#B2_FILE="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/txt/Macrophage_tracks.txt"

GTF_FILE="/bobross/STAR_Metadata/mm10/gencode.vM23.annotation.gtf"
READ_LENGTH=38
NTHREAD=16
OUTPUT_DIR="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_only"
#OUTPUT_DIR="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/T_only"
#OUTPUT_DIR="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/Macrophage_only"

TMP_DIR="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/tmp"
SEQ_TYPE="paired"

# Set the library path
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# Run rmats.py with the specified parameters
python "$RMATS_FILE" \
    --b1 "$B1_FILE" \
    --b2 "$B2_FILE" \
    --gtf "$GTF_FILE" \
    -t "$SEQ_TYPE" \
    --readLength "$READ_LENGTH" \
    --nthread "$NTHREAD" \
    --od "$OUTPUT_DIR" \
    --tmp "$TMP_DIR"
