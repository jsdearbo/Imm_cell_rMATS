#!/usr/bin/env bash
set -euo pipefail

########################
# USER SETTINGS
########################

GUIDE_GTF="/bobross/STAR_Metadata/mm10/gencode.vM23.annotation.gtf"
BAM_DIR="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/bam"
OUTDIR="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/compute_psi/stringtie_out"
THREADS=8

# strandedness
#STRAND_OPT="--rf"   # fr-firststrand
#STRAND_OPT="--fr"   # fr-secondstrand
STRAND_OPT=""        # unstranded

########################
# 1) PER-SAMPLE ASSEMBLY
########################

mkdir -p "${OUTDIR}"

echo "[INFO] Running StringTie per sample..."

# ensure unmatched globs just disappear instead of literally 'B_Cell*.bam'
shopt -s nullglob

# ONLY run on these three groups
for BAM in "${BAM_DIR}"/B_Cell*.bam \
           "${BAM_DIR}"/T_Cell*.bam \
           "${BAM_DIR}"/Macrophage*.bam; do

    SAMPLE=$(basename "${BAM}" .bam)
    echo "[INFO] Processing ${SAMPLE}"

    stringtie "${BAM}" \
        -G "${GUIDE_GTF}" \
        -o "${OUTDIR}/${SAMPLE}.gtf" \
        -p "${THREADS}" \
        ${STRAND_OPT} \
        -l "${SAMPLE}"
done

########################
# 2) CELL-TYPE–SPECIFIC MERGES
########################

echo "[INFO] Creating cell-type–specific GTF lists..."

B_PREFIX="B_Cell"
T_PREFIX="T_Cell"
M_PREFIX="Macrophage"

B_LIST="${OUTDIR}/B_Cell_gtf_list.txt"
T_LIST="${OUTDIR}/T_Cell_gtf_list.txt"
M_LIST="${OUTDIR}/Macrophage_gtf_list.txt"

ls "${OUTDIR}"/${B_PREFIX}*.gtf  > "${B_LIST}"
ls "${OUTDIR}"/${T_PREFIX}*.gtf  > "${T_LIST}"
ls "${OUTDIR}"/${M_PREFIX}*.gtf  > "${M_LIST}"

echo "[INFO] Merging B cell transcripts..."
stringtie --merge \
    -G "${GUIDE_GTF}" \
    -o "${OUTDIR}/B_Cell_merged.gtf" \
    "${B_LIST}"

echo "[INFO] Merging T cell transcripts..."
stringtie --merge \
    -G "${GUIDE_GTF}" \
    -o "${OUTDIR}/T_Cell_merged.gtf" \
    "${T_LIST}"

echo "[INFO] Merging Macrophage transcripts..."
stringtie --merge \
    -G "${GUIDE_GTF}" \
    -o "${OUTDIR}/Macrophage_merged.gtf" \
    "${M_LIST}"

########################
# 3) GLOBAL MERGE (ONLY THESE 3 TYPES)
########################

ALL_LIST="${OUTDIR}/all_gtf_list.txt"

# only include B_Cell, T_Cell, Macrophage; avoid any old gd_T_Cell*.gtf etc.
ls "${OUTDIR}"/B_Cell*.gtf \
   "${OUTDIR}"/T_Cell*.gtf \
   "${OUTDIR}"/Macrophage*.gtf \
   > "${ALL_LIST}"

echo "[INFO] Merging ALL samples (B/T/M only)..."
stringtie --merge \
    -G "${GUIDE_GTF}" \
    -o "${OUTDIR}/ALL_3_merged.gtf" \
    "${ALL_LIST}"

echo "[INFO] Done."
echo "[INFO] B-cell GTF:        ${OUTDIR}/B_Cell_merged.gtf"
echo "[INFO] T-cell GTF:        ${OUTDIR}/T_Cell_merged.gtf"
echo "[INFO] Macrophage GTF:    ${OUTDIR}/Macrophage_merged.gtf"
echo "[INFO] Global union GTF:  ${OUTDIR}/ALL_3_merged.gtf"
