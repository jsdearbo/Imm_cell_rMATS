#!/usr/bin/env bash
set -euo pipefail

# Paths to text files produced for rMATS input (comma- or newline-delimited BAM paths)
B_CELL_BAM_LIST=/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/B_Cell_bams.txt
T_CELL_BAM_LIST=/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/T_Cell_bams.txt
MACROPHAGE_BAM_LIST=/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/Macrophage_bams.txt

# Output BAM destinations
B_CELL_MERGED_BAM=/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/merge_bams/B_Cell_merged.bam
T_CELL_MERGED_BAM=/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/merge_bams/T_Cell_merged.bam
MACROPHAGE_MERGED_BAM=/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/merge_bams/Macrophage_merged.bam

THREADS=${THREADS:-4}

require_tools() {
    command -v samtools >/dev/null 2>&1 || {
        echo "samtools not on PATH; activate the environment first" >&2
        exit 127
    }
}

merge_bams() {
    local bam_list_file=$1
    local output_bam=$2

    mkdir -p "$(dirname "$output_bam")"

    local tmp_list
    tmp_list=$(mktemp)
    trap 'rm -f "${tmp_list:-}"' RETURN

    tr ',' '\n' < "$bam_list_file" \
        | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' \
        | sed '/^$/d' \
        > "$tmp_list"

    if [[ ! -s $tmp_list ]]; then
        echo "No BAM paths found in $bam_list_file" >&2
        return 1
    fi

    # Merge and preserve headers/program groups from all inputs
    samtools merge -@ "$THREADS" -c -p -b "$tmp_list" "$output_bam"

    samtools index -@ "$THREADS" "$output_bam"
}

main() {
    require_tools

    merge_bams "$B_CELL_BAM_LIST" "$B_CELL_MERGED_BAM"
    merge_bams "$T_CELL_BAM_LIST" "$T_CELL_MERGED_BAM"
    merge_bams "$MACROPHAGE_BAM_LIST" "$MACROPHAGE_MERGED_BAM"
}

main "$@"