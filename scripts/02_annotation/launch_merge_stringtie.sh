#!/usr/bin/env bash
set -euo pipefail

# Run BAM merge workflow followed by StringTie annotation generation.

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
merge_script="${script_dir}/merge_cell_type_bams.sh"
stringtie_script="${script_dir}/stringtie_annotation.sh"

if [[ ! -x "$merge_script" ]]; then
	echo "Missing or non-executable merge script: $merge_script" >&2
	exit 1
fi

if [[ ! -x "$stringtie_script" ]]; then
	echo "Missing or non-executable StringTie script: $stringtie_script" >&2
	exit 1
fi

# Reproducible StringTie defaults (can still be overridden by caller env)
: "${STRICT_REFERENCE:=1}"
: "${MIN_COVERAGE:=2.0}"
: "${MIN_JUNCTION_COV:=2}"
: "${MIN_ISOFORM_FRAC:=0.05}"
: "${MIN_TRANSCRIPT_LEN:=200}"

export STRICT_REFERENCE
export MIN_COVERAGE
export MIN_JUNCTION_COV
export MIN_ISOFORM_FRAC
export MIN_TRANSCRIPT_LEN

echo "[INFO] StringTie settings: STRICT_REFERENCE=${STRICT_REFERENCE} MIN_COVERAGE=${MIN_COVERAGE} MIN_JUNCTION_COV=${MIN_JUNCTION_COV} MIN_ISOFORM_FRAC=${MIN_ISOFORM_FRAC} MIN_TRANSCRIPT_LEN=${MIN_TRANSCRIPT_LEN}"

"$merge_script"

"$stringtie_script"
