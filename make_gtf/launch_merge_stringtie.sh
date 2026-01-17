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

"$merge_script"

"$stringtie_script"
