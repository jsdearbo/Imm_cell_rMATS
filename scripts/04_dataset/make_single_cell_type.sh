#!/usr/bin/env bash
# ==============================================================================
# make_single_cell_type.sh - Extract single-cell-type NPZ datasets from the
# multitask NPZ shards produced by Step 3 (launch_convert_parallel.sh).
#
# For each cell type, filters targets/mask/confidence to a single task column,
# producing train/val/test NPZ shards in SINGLE_CELL_OUT_DIR/{cell_type}/.
#
# Run AFTER Step 3 is complete and multitask_balanced shards exist.
#
# Usage:
#   bash scripts/04_dataset/make_single_cell_type.sh
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG="${REPO_ROOT}/config/pipeline_config.sh"

if [[ ! -f "${CONFIG}" ]]; then
    echo "[ERROR] Missing config: ${CONFIG}" >&2
    exit 1
fi
source "${CONFIG}"

echo "============================================================"
echo "  [Step 4] Single-cell NPZ extraction"
echo "  Input:  ${NPZ_OUT_DIR}"
echo "  Output: ${SINGLE_CELL_OUT_DIR}"
echo "============================================================"

for CELL_TYPE in "${CELL_TYPES_ARR[@]}"; do
    OUT_DIR="${SINGLE_CELL_OUT_DIR}/${CELL_TYPE}"
    mkdir -p "${OUT_DIR}"
    echo "[Step 4] Processing cell type: ${CELL_TYPE} -> ${OUT_DIR}"
    python "${FILTER_SCRIPT}" \
        --input-dir  "${NPZ_OUT_DIR}" \
        --cell-type  "${CELL_TYPE}" \
        --output-dir "${OUT_DIR}"
done

echo ""
echo "============================================================"
echo "  [Step 4] Done.  Single-cell datasets in:"
echo "  ${SINGLE_CELL_OUT_DIR}/"
echo "============================================================"
