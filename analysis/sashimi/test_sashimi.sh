#!/usr/bin/env bash
set -euo pipefail

############################
# USER CONFIG
############################

# If rmats2sashimiplot is installed as a console script, this can just be "rmats2sashimiplot"
RMATS2SASHIMI_PY="/mnt/netfiles/jsdearbo/miniconda3/envs/rmats/bin/rmats2sashimiplot"

# rMATS output directory (contains SE.MATS.JC.txt, etc.)
RMATS_OUT_DIR="/mnt/netfiles/jsdearbo/public_data/ImmGen_Deep_RNAseq/rmats/B_vs_mac"

# Plot labels
LABEL1="B cell"
LABEL2="Macrophage"

# BAM lists used in rMATS (newline-separated)
B1_LIST="//mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/B_Cell_tracks.txt"
B2_LIST="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/txt/Macrophage_tracks.txt"

# Output
SASHIMI_OUT_BASE="/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/B_vs_mac/sashimi_B_vs_mac"
TMPDIR="${SASHIMI_OUT_BASE}/tmp_lists"
mkdir -p "$TMPDIR"

# chr-prefix handling for events vs BAMs: keep | remove | auto
CHR_PREFIX_MODE="auto"

# Visual tweaks
EXON_SCALE=1
INTRON_SCALE=8
MIN_COUNTS=5
FONT_SIZE=9
FIG_HEIGHT=""   # e.g., 7 to force; leave blank for default
FIG_WIDTH=""    # e.g., 10 to force; leave blank for default

# Filtering
FDR_THRESH=0.05          # e.g., 0.05 for significant-only; 1.0 keeps all tested events
ABS_DPSI_MIN=0.5          # e.g., 0.1 to require |ΔPSI| >= 0.1
MAX_EVENTS=500          # cap to avoid marathon runs

# Which event types to run
event_types=(SE A5SS A3SS MXE)   # change to: (SE A5SS A3SS MXE RI)

############################
# PREP
############################

# Build comma-separated BAM lists
b1_csv="${TMPDIR}/b1.csv"
b2_csv="${TMPDIR}/b2.csv"
paste -sd, <(sed '/^\s*$/d;s/^[ \t]*//;s/[ \t]*$//' "$B1_LIST") > "$b1_csv"
paste -sd, <(sed '/^\s*$/d;s/^[ \t]*//;s/[ \t]*$//' "$B2_LIST") > "$b2_csv"

# Count replicates
n1=$(( $(tr -cd ',' < "$b1_csv" | wc -c) + 1 ))
n2=$(( $(tr -cd ',' < "$b2_csv" | wc -c) + 1 ))

# Build VALID group file using 1-based indices (over concatenated b1 then b2)
gf="${TMPDIR}/groups.gf"
start2=$((n1+1))
end2=$((n1+n2))
printf "Bcell: 1-%d\nMac: %d-%d\n" "$n1" "$start2" "$end2" > "$gf"
echo "[OK] Group file:"
cat "$gf"

# Optional args
opt_args=()
[[ -n "${EXON_SCALE}"  ]] && opt_args+=( --exon_s "${EXON_SCALE}" )
[[ -n "${INTRON_SCALE}" ]] && opt_args+=( --intron_s "${INTRON_SCALE}" )
[[ -n "${MIN_COUNTS}"  ]] && opt_args+=( --min-counts "${MIN_COUNTS}" )
[[ -n "${FONT_SIZE}"   ]] && opt_args+=( --font-size "${FONT_SIZE}" )
[[ -n "${FIG_HEIGHT}"  ]] && opt_args+=( --fig-height "${FIG_HEIGHT}" )
[[ -n "${FIG_WIDTH}"   ]] && opt_args+=( --fig-width  "${FIG_WIDTH}" )

case "$CHR_PREFIX_MODE" in
  keep)   opt_args+=( --keep-event-chr-prefix ) ;;
  remove) opt_args+=( --remove-event-chr-prefix ) ;;
  auto)   : ;;  # no flag
  *) echo "Invalid CHR_PREFIX_MODE: $CHR_PREFIX_MODE" >&2; exit 1 ;;
esac

############################
# RUN
############################

echo "[INFO] Launching rmats2sashimiplot for: ${event_types[*]}"
for et in "${event_types[@]}"; do
  events_file="${RMATS_OUT_DIR}/${et}.MATS.JC.txt"

  if [[ ! -s "$events_file" ]]; then
    echo "[WARN] Missing rMATS event file for ${et}: ${events_file} — skipping."
    continue
  fi

  out_dir="${SASHIMI_OUT_BASE}/${et}"
  mkdir -p "$out_dir"

  # -------- Filter rMATS events before plotting --------
  subset="${TMPDIR}/${et}.filtered.txt"

  # Count original events (excluding header) for reporting
  orig_count=$(( $(wc -l < "$events_file") - 1 ))
  (( orig_count < 0 )) && orig_count=0

  awk -v fdr="$FDR_THRESH" -v dpsi="$ABS_DPSI_MIN" -v maxn="$MAX_EVENTS" '
    BEGIN {
      FS=OFS="\t"
      n=0
      BIG=1e9
    }
    NR==1 {
      # locate columns by name (robust to column order)
      for (i=1;i<=NF;i++){
        if ($i=="FDR") fdr_col=i
        if ($i=="IncLevelDifference") dpsi_col=i
        if ($i=="IncLevel1") il1_col=i
        if ($i=="IncLevel2") il2_col=i
        if ($i=="chr") chr_col=i
      }
      print; next
    }
    {
      # chromosome (skip sex and mitochondrial)
      c = (chr_col ? $chr_col : $1)
      if (c ~ /^chr[XY]$/) next
      if (c == "chrM" || c == "MT") next

      # skip if any IncLevel column contains NA
      if (il1_col && $il1_col ~ /NA/) next
      if (il2_col && $il2_col ~ /NA/) next

      # numeric guards
      f  = (fdr_col  ? $fdr_col  + 0 : BIG)
      dp = (dpsi_col ? $dpsi_col + 0 : 0)

      # FDR and |ΔPSI| filters
      if (f < fdr && (dpsi==0 || (dp >= dpsi || dp <= -dpsi))) {
        print
        n++
        if (n>=maxn) exit
      }
    }
  ' "$events_file" > "$subset"

  # Count remaining (excluding header) and report
  rem_count=$(( $(wc -l < "$subset") - 1 ))
  (( rem_count < 0 )) && rem_count=0
  echo "[INFO] ${et}: ${rem_count}/${orig_count} events after filters (FDR<${FDR_THRESH}, |ΔPSI|>=${ABS_DPSI_MIN}, no chrX/Y/M)."

  # If nothing survived filtering, skip plotting this event type
  if [[ "$rem_count" -le 0 ]]; then
    echo "[INFO] ${et}: no events passed filters. Skipping."
    continue
  fi

  events_file="$subset"
  # ------------------------------------------------------

  echo "[RUN] ${et} -> ${out_dir}"
  "$RMATS2SASHIMI_PY" \
    -o "${out_dir}" \
    --l1 "${LABEL1}" \
    --l2 "${LABEL2}" \
    --event-type "${et}" \
    -e "${events_file}" \
    --b1 "${b1_csv}" \
    --b2 "${b2_csv}" \
    --group-info "${gf}" \
    "${opt_args[@]}"
done

echo "[DONE] Sashimi plots written under: ${SASHIMI_OUT_BASE}"
