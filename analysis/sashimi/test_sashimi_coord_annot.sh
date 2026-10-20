#!/usr/bin/env bash
set -euo pipefail

############################
# USER CONFIG
############################

RMATS2SASHIMI="/bobross/miniconda3/envs/rmats/bin/rmats2sashimiplot"
RMATS_OUT_DIR="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_vs_mac"
GFF3="/bobross/STAR_Metadata/mm10/gencode.vM23.annotation.gff3"

LABEL1="B cell"
LABEL2="Macrophage"

B1_LIST="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/txt/B_Cell_tracks.txt"
B2_LIST="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/txt/Macrophage_tracks.txt"

OUT_BASE="/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/sashimi_coord_mode"
TMPDIR="${OUT_BASE}/tmp"; mkdir -p "$TMPDIR"

# Visual tweaks
EXON_SCALE=1
INTRON_SCALE=8
MIN_COUNTS=5
FONT_SIZE=9

# Filtering
FDR_THRESH=0.5
ABS_DPSI_MIN=0
MAX_EVENTS=100
PAD_BP=200

# Event types to run
event_types=(SE)

############################
# PREP
############################

b1_csv="${TMPDIR}/b1.csv"
b2_csv="${TMPDIR}/b2.csv"
paste -sd, <(sed '/^\s*$/d;s/^[ \t]*//;s/[ \t]*$//' "$B1_LIST") > "$b1_csv"
paste -sd, <(sed '/^\s*$/d;s/^[ \t]*//;s/[ \t]*$//' "$B2_LIST") > "$b2_csv"

n1=$(( $(tr -cd ',' < "$b1_csv" | wc -c) + 1 ))
n2=$(( $(tr -cd ',' < "$b2_csv" | wc -c) + 1 ))
gf="${TMPDIR}/groups.gf"
printf "Bcell: 1-%d\nMac: %d-%d\n" "$n1" "$((n1+1))" "$((n1+n2))" > "$gf"
echo "[OK] Group file:"; cat "$gf"

plot_args=( --exon_s "$EXON_SCALE" --intron_s "$INTRON_SCALE" --min-counts "$MIN_COUNTS" --font-size "$FONT_SIZE" )

############################
# RUN (coordinate + annotation)
############################

for et in "${event_types[@]}"; do
  events="${RMATS_OUT_DIR}/${et}.MATS.JC.txt"
  [[ -s "$events" ]] || { echo "[WARN] Missing $events, skip."; continue; }

  out_dir="${OUT_BASE}/${et}"; mkdir -p "$out_dir"

  # Filter rMATS table (robust to column order)
  filtered="${TMPDIR}/${et}.filtered.txt"
  orig_count=$(( $(wc -l < "$events") - 1 )); (( orig_count < 0 )) && orig_count=0

  awk -v fdr="$FDR_THRESH" -v dpsi="$ABS_DPSI_MIN" -v maxn="$MAX_EVENTS" '
    BEGIN { FS=OFS="\t"; BIG=1e9; n=0 }
    NR==1 {
      for (i=1;i<=NF;i++){
        key=$i
        gsub(/\r$/,"",key)         # strip CR if present
        if (key=="chr") chr_col=i
        if (key=="strand") strand_col=i
        if (key=="upstreamES") uES=i
        if (key=="upstreamEE") uEE=i
        if (key=="exonStart_0base") exS=i
        if (key=="exonEnd") exE=i
        if (key=="downstreamES") dES=i
        if (key=="downstreamEE") dEE=i
        if (key=="FDR") fdr_col=i
        if (key=="IncLevelDifference") dpsi_col=i
        if (key=="IncLevel1") il1_col=i
        if (key=="IncLevel2") il2_col=i
      }
      print; next
    }
    {
      c  = (chr_col     ? $chr_col     : $1)
      st = (strand_col  ? $strand_col  : "+")
      if (c ~ /^chr[XY]$/) next
      if (c == "chrM" || c == "MT") next
      if (il1_col && $il1_col ~ /NA/) next
      if (il2_col && $il2_col ~ /NA/) next

      F  = (fdr_col  ? $fdr_col  + 0 : BIG)
      DP = (dpsi_col ? $dpsi_col + 0 : 0)

      if (F < fdr && (dpsi==0 || (DP >= dpsi || DP <= -dpsi))) {
        print
        n++; if (n>=maxn) exit
      }
    }
  ' "$events" > "$filtered"

  rem_count=$(( $(wc -l < "$filtered") - 1 )); (( rem_count < 0 )) && rem_count=0
  echo "[INFO] ${et}: ${rem_count}/${orig_count} events after filters (FDR<${FDR_THRESH}, |ΔPSI|>=${ABS_DPSI_MIN}, no chrX/Y/M)."
  [[ "$rem_count" -gt 0 ]] || { echo "[INFO] ${et}: nothing to plot after filters."; continue; }

  # Build robust coords for SE
  coords="${TMPDIR}/${et}.coords.txt"
  awk -v PAD="$PAD_BP" '
    BEGIN { FS=OFS="\t" }
    NR==1 {
      for (i=1;i<=NF;i++){
        key=$i; gsub(/\r$/,"",key)
        if (key=="chr") chr_col=i
        if (key=="strand") strand_col=i
        if (key=="upstreamES") uES=i
        if (key=="upstreamEE") uEE=i
        if (key=="exonStart_0base") exS=i
        if (key=="exonEnd") exE=i
        if (key=="downstreamES") dES=i
        if (key=="downstreamEE") dEE=i
      }
      next
    }
    {
      c  = (chr_col    ? $chr_col    : $1)
      st = (strand_col ? $strand_col : "+")
      # prefer window spanning flanking exons; if any missing, fall back to middle exon
      s = (uES ? $uES + 0 : (exS ? $exS + 0 : 0))
      e = (dEE ? $dEE + 0 : (exE ? $exE + 0 : 0))

      # ensure s<e; if reversed for any reason, swap
      if (s>e) { t=s; s=e; e=t }

      if (s>0 && e>0 && e>s) {
        s = (s - PAD < 0 ? 0 : s - PAD)
        e = e + PAD
        print c, st, s, e
      }
    }
  ' "$filtered" > "$coords"

  # Sanity: show a preview
  echo "[DEBUG] First 5 coord lines:"
  head -n 5 "$coords" || true

  # Quick chrom-name check against BAM headers (first BAM only)
  first_bam=$(cut -d, -f1 "$b1_csv")
  if command -v samtools >/dev/null 2>&1; then
    first_chrom_in_bam=$(samtools idxstats "$first_bam" | awk 'NR==1{print $1}')
    first_chrom_in_gff=$(awk '$1!~/#/ {print $1; exit}' "$GFF3")
    echo "[DEBUG] BAM first chrom: $first_chrom_in_bam ; GFF3 first chrom: $first_chrom_in_gff"
  else
    echo "[DEBUG] samtools not found; skipping chrom-name cross-check."
  fi

  # Guard: no coords?
  if [[ ! -s "$coords" ]]; then
    echo "[WARN] No coordinates derived for ${et}. Check header names in ${filtered}."
    continue
  fi

  # Iterate coords and plot; don’t die on a single failure
  idx=0
  while IFS=$'\t' read -r chrom strand start end; do
    ((idx++))
    region="${chrom}:${strand}:${start}:${end}:${GFF3}"
    ev_out="${out_dir}/${chrom}_${start}_${end}_${strand}"
    mkdir -p "$ev_out"

    echo "[RUN] ${et} #${idx}: ${region}"
    set +e
    "$RMATS2SASHIMI" \
      -o "$ev_out" \
      --l1 "$LABEL1" \
      --l2 "$LABEL2" \
      -c "$region" \
      --b1 "$b1_csv" \
      --b2 "$b2_csv" \
      --group-info "$gf" \
      "${plot_args[@]}"
    rc=$?
    set -e
    if [[ $rc -ne 0 ]]; then
      echo "[WARN] Plot failed for ${region} (exit ${rc}); continuing."
      continue
    fi
  done < "$coords"

done

echo "[DONE] Coordinate+annotation sashimi written under: ${OUT_BASE}"
