import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import pyranges as pr

from io_utils import load_label_intervals
from generate_splicing_labels import (
    load_rmats_se,
    load_rmats_ri,
    load_rmats_alt_site,
    load_stringtie_gtf,
    generate_soft_labels,
    collapse_exon_coords_weighted,
    _disjoint_labels_with_priority,
)

# ---- Project layout ----
_PROJECT_ROOT = "/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS"
_RESULTS_ROOT = f"{_PROJECT_ROOT}/results"
_LABELING_ROOT = f"{_PROJECT_ROOT}/data_labeling"

DEFAULTS = {
    # Single-mode defaults (B_Cell as example)
    "rmats_se":    f"{_RESULTS_ROOT}/B_only/SE.MATS.JC.txt",
    "rmats_ri":    f"{_RESULTS_ROOT}/B_only/RI.MATS.JC.txt",
    "rmats_a3ss":  f"{_RESULTS_ROOT}/B_only/A3SS.MATS.JC.txt",
    "rmats_a5ss":  f"{_RESULTS_ROOT}/B_only/A5SS.MATS.JC.txt",
    "gtf":         f"{_RESULTS_ROOT}/stringtie_per_rep/dominant_isoforms/B_Cell_dominant.gtf",
    "min_coverage": 20,
    "out_bed":     f"{_LABELING_ROOT}/test_soft_labeling_output/B_cell_soft_labels.bed",
    "summary":     f"{_LABELING_ROOT}/test_soft_labeling_output/B_cell_soft_label_summary.tsv",
    # Master-mode defaults
    "rmats_root":        _RESULTS_ROOT,
    "gtf_root":          f"{_RESULTS_ROOT}/stringtie_per_rep/dominant_isoforms",
    "per_sample_out_dir": f"{_LABELING_ROOT}/outputs_soft_per_sample/soft",
    "cell_types":        "B_Cell,T_Cell,Macrophage",
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate soft-label BEDs from rMATS output, with optional master BED creation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Two modes of operation:

1. SINGLE MODE (default): Generate one BED file from rMATS files
   python test_soft_labeling.py --rmats-se SE.txt --rmats-ri RI.txt --gtf dominant.gtf --out-bed output.bed

2. MASTER MODE: Create master BED files by aggregating replicates per cell type
   python test_soft_labeling.py --manifest manifest.tsv --out-dir /path/to/output

   Manifest format (tab-separated):
   cell_type    bed_path
   B_cell       /path/to/B_cell_rep1.bed
   B_cell       /path/to/B_cell_rep2.bed
   T_cell       /path/to/T_cell_rep1.bed
        """
    )
    
    # Mode selection
    parser.add_argument("--manifest", help="Optional TSV manifest for master BED mode (columns must include cell_type, bed_path)")
    parser.add_argument("--out-dir", help="Output directory for master BED mode (master_beds/ will be created inside)")
    parser.add_argument("--regenerate-beds", dest="regenerate_beds", action="store_true", help="In manifest mode, regenerate per-sample BEDs from rMATS before master aggregation (default)")
    parser.add_argument("--no-regenerate-beds", dest="regenerate_beds", action="store_false", help="In manifest mode, skip regeneration and aggregate existing BEDs only")
    parser.set_defaults(regenerate_beds=True)
    parser.add_argument("--rmats-root", default=DEFAULTS["rmats_root"], help="Root directory containing <celltype>_only rMATS subdirs")
    parser.add_argument("--gtf-root", default=DEFAULTS["gtf_root"], help="Directory containing <celltype>_dominant.gtf files")
    parser.add_argument("--per-sample-out-dir", default=DEFAULTS["per_sample_out_dir"], help="Output directory for regenerated per-sample BEDs")
    parser.add_argument("--manifest-out", default=None, help="Optional path to write an updated manifest after regeneration")
    parser.add_argument("--cell-types", default=DEFAULTS["cell_types"], help="Comma-separated cell types for auto-manifest mode (used when --manifest is omitted)")
    
    # Single mode arguments (used when --manifest not provided)
    parser.add_argument("--rmats-se", default=DEFAULTS["rmats_se"], help="Path to SE.MATS.JC.txt")
    parser.add_argument("--rmats-ri", default=DEFAULTS["rmats_ri"], help="Path to RI.MATS.JC.txt")
    parser.add_argument("--rmats-a3ss", default=DEFAULTS["rmats_a3ss"], help="Path to A3SS.MATS.JC.txt")
    parser.add_argument("--rmats-a5ss", default=DEFAULTS["rmats_a5ss"], help="Path to A5SS.MATS.JC.txt")
    parser.add_argument("--gtf", default=DEFAULTS["gtf"], help="Dominant StringTie GTF")
    parser.add_argument("--min-coverage", type=int, default=DEFAULTS["min_coverage"], help="rMATS coverage cutoff")
    parser.add_argument("--out-bed", default=DEFAULTS["out_bed"], help="Output BED file (chr start end p_exon confidence)")
    parser.add_argument("--summary", default=DEFAULTS["summary"], help="Optional TSV summary of counts by label")
    
    return parser.parse_args()

def create_master_bed_files(bed_files_by_celltype: dict, out_dir: str):
    """
    Creates aggregated master BED files for each cell type.
        Aggregates replicates for matching intervals using read-depth weighting.

        Weight source:
            confidence ~= log1p(read_depth), so read_depth ~= expm1(confidence)

        For each (Chromosome, Start, End):
            PSI_master = sum(PSI_i * depth_i) / sum(depth_i)
            confidence_master = log1p(sum(depth_i))

        If all depths are zero/missing for an interval, falls back to unweighted PSI mean.
    """
    master_dir = Path(out_dir) / "master_beds"
    master_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"[INFO] Creating master BED files in {master_dir}")
    
    for ct, bed_paths in bed_files_by_celltype.items():
        print(f"[INFO]   Processing {ct}: {len(bed_paths)} replicate(s)")
        
        # Load all intervals from all replicates
        all_intervals = []
        for bpath in bed_paths:
            df = load_label_intervals(bpath)
            if df.empty:
                print(f"[WARN]     Skipping empty BED file: {bpath}")
                continue
            all_intervals.append(df)
        
        if not all_intervals:
            print(f"[WARN]   No valid intervals for {ct}, skipping.")
            continue
        
        # Concatenate all intervals
        combined = pd.concat(all_intervals, ignore_index=True)
        print(f"[INFO]     Loaded {len(combined)} total intervals from all replicates")
        
        # Sort by chromosome and position
        combined = combined.sort_values(["Chromosome", "Start", "End"]).reset_index(drop=True)

        # Derive read-depth proxy from confidence (confidence = log1p(depth))
        combined["depth"] = np.expm1(combined["confidence"].clip(lower=0.0)).astype(float)
        combined["weighted_psi"] = combined["PSI"].astype(float) * combined["depth"]

        grouped = combined.groupby(["Chromosome", "Start", "End"], as_index=False).agg(
            psi_num=("weighted_psi", "sum"),
            depth_sum=("depth", "sum"),
            psi_mean=("PSI", "mean"),
        )

        grouped["PSI"] = np.where(
            grouped["depth_sum"] > 0,
            grouped["psi_num"] / grouped["depth_sum"],
            grouped["psi_mean"],
        )
        grouped["confidence"] = np.log1p(grouped["depth_sum"]).astype(float)
        grouped = grouped[["Chromosome", "Start", "End", "PSI", "confidence"]]
        
        # Write master BED file
        output_path = master_dir / f"{ct}_master.bed"
        grouped.to_csv(
            output_path,
            sep="\t",
            header=False,
            index=False,
            columns=["Chromosome", "Start", "End", "PSI", "confidence"]
        )
        
        print(f"[INFO]     Wrote {len(grouped)} unique intervals to {output_path.name}")
    
    print("[INFO] Master BED file generation complete.")


def get_psi_sample_columns(df: pd.DataFrame):
    psi_cols = []
    for c in df.columns:
        if not c.startswith("psi_"):
            continue
        suffix = c.split("_", 1)[1]
        if suffix.isdigit():
            psi_cols.append(c)
    return sorted(psi_cols, key=lambda x: int(x.split("_")[1]))


def _celltype_rmats_prefix(cell_type: str) -> str:
    if cell_type in {"B_Cell", "T_Cell"}:
        return cell_type.split("_")[0]
    return cell_type


def _parse_cell_types(cell_types_arg: str):
    return [c.strip() for c in str(cell_types_arg).split(",") if c.strip()]


def _resolve_celltype_paths(cell_type: str, args, manifest_sub: pd.DataFrame = None):
    rmats_prefix = _celltype_rmats_prefix(cell_type)
    if manifest_sub is None:
        manifest_sub = pd.DataFrame()

    rmats_se = manifest_sub["rmats_se"].iloc[0] if "rmats_se" in manifest_sub.columns and pd.notna(manifest_sub["rmats_se"].iloc[0]) else str(Path(args.rmats_root) / f"{rmats_prefix}_only" / "SE.MATS.JC.txt")
    rmats_ri = manifest_sub["rmats_ri"].iloc[0] if "rmats_ri" in manifest_sub.columns and pd.notna(manifest_sub["rmats_ri"].iloc[0]) else str(Path(args.rmats_root) / f"{rmats_prefix}_only" / "RI.MATS.JC.txt")
    rmats_a3ss = manifest_sub["rmats_a3ss"].iloc[0] if "rmats_a3ss" in manifest_sub.columns and pd.notna(manifest_sub["rmats_a3ss"].iloc[0]) else str(Path(args.rmats_root) / f"{rmats_prefix}_only" / "A3SS.MATS.JC.txt")
    rmats_a5ss = manifest_sub["rmats_a5ss"].iloc[0] if "rmats_a5ss" in manifest_sub.columns and pd.notna(manifest_sub["rmats_a5ss"].iloc[0]) else str(Path(args.rmats_root) / f"{rmats_prefix}_only" / "A5SS.MATS.JC.txt")
    gtf = manifest_sub["gtf"].iloc[0] if "gtf" in manifest_sub.columns and pd.notna(manifest_sub["gtf"].iloc[0]) else str(Path(args.gtf_root) / f"{cell_type}_dominant.gtf")

    return rmats_se, rmats_ri, rmats_a3ss, rmats_a5ss, gtf


def regenerate_one_cell_type(cell_type: str, args, sub: pd.DataFrame = None):
    if sub is None:
        sub = pd.DataFrame()

    print(f"[INFO] Regenerating per-sample BEDs for {cell_type}")

    rmats_se_p, rmats_ri_p, rmats_a3ss_p, rmats_a5ss_p, gtf_p = _resolve_celltype_paths(cell_type, args, sub)
    print(f"[INFO]   rMATS dir inputs: {rmats_se_p}, {rmats_ri_p}, {rmats_a3ss_p}, {rmats_a5ss_p}")
    print(f"[INFO]   GTF: {gtf_p}")

    if not Path(rmats_se_p).exists():
        print(f"[WARN]   Missing SE file for {cell_type}: {rmats_se_p}; skipping")
        return [], []
    if not Path(gtf_p).exists():
        print(f"[WARN]   Missing GTF for {cell_type}: {gtf_p}; skipping")
        return [], []

    rmats_se = load_rmats_se(rmats_se_p, args.min_coverage)
    rmats_ri = load_rmats_ri(rmats_ri_p, args.min_coverage)
    rmats_a3ss = load_rmats_alt_site(rmats_a3ss_p, args.min_coverage)
    rmats_a5ss = load_rmats_alt_site(rmats_a5ss_p, args.min_coverage)
    stringtie_gr = load_stringtie_gtf(gtf_p)

    vars_list = [gr for gr in (rmats_se, rmats_ri, rmats_a3ss, rmats_a5ss) if len(gr) > 0]
    if not vars_list:
        print(f"[WARN]   No rMATS events loaded for {cell_type}; skipping regeneration")
        return [], []

    rmats_all = pr.concat(vars_list)
    if len(rmats_all) == 0:
        print(f"[WARN]   Empty concatenated rMATS intervals for {cell_type}; skipping")
        return [], []

    segments = build_shared_segments_df(stringtie_gr, rmats_all.df)
    collapsed = collapse_exon_coords_weighted(rmats_all.df)
    psi_cols = get_psi_sample_columns(collapsed)
    if not psi_cols:
        print(f"[WARN]   No per-sample PSI columns after collapse for {cell_type}; skipping")
        return [], []

    psi_to_out = _map_psi_to_output_paths(cell_type, sub, psi_cols, args.per_sample_out_dir)

    produced_rows = []
    for psi_col in psi_cols:
        psi_events = (
            collapsed[["Chromosome", "Start", "End", "Strand", psi_col]]
            .rename(columns={psi_col: "psi"})
            .dropna(subset=["psi"])
        )

        bed = _disjoint_labels_with_priority(
            psi_df=psi_events,
            segments_df=segments,
        )

        sample_idx = int(psi_col.split("_")[1])
        ijc_col = f"ijc_{sample_idx}"
        sjc_col = f"sjc_{sample_idx}"
        if ijc_col in collapsed.columns and sjc_col in collapsed.columns:
            conf_events = collapsed[["Chromosome", "Start", "End", "Strand", ijc_col, sjc_col]].copy()
            conf_events["confidence"] = np.log1p(conf_events[ijc_col].astype(float) + conf_events[sjc_col].astype(float))
            conf_events = conf_events[["Chromosome", "Start", "End", "Strand", "confidence"]]
        else:
            conf_events = pd.DataFrame(columns=["Chromosome", "Start", "End", "Strand", "confidence"])

        conf_on_seg = _project_values_to_segments(segments, conf_events, "confidence", default_value=0.0)
        bed = bed.merge(conf_on_seg, on=["Chromosome", "Start", "End", "Strand"], how="left")
        bed["confidence"] = bed["confidence"].fillna(0.0)

        out_path = Path(psi_to_out[psi_col])
        out_path.parent.mkdir(parents=True, exist_ok=True)
        bed[["Chromosome", "Start", "End", "p_exon", "confidence"]].to_csv(
            out_path,
            sep="\t",
            header=False,
            index=False,
        )
        print(f"[INFO]   Wrote {out_path}")

        produced_rows.append({
            "sample_id": f"{cell_type}_{psi_col}",
            "cell_type": cell_type,
            "bed_path": str(out_path),
            "replicate_id": str(sample_idx),
            "rmats_se": rmats_se_p,
            "rmats_ri": rmats_ri_p,
            "rmats_a3ss": rmats_a3ss_p,
            "rmats_a5ss": rmats_a5ss_p,
            "gtf": gtf_p,
        })

    return produced_rows, psi_cols


def build_shared_segments_df(stringtie_gr, rmats_all_bounds_df):
    gdf = stringtie_gr.df
    exons_df = gdf[gdf.Feature == "exon"][["Chromosome", "Start", "End", "Strand"]].copy()
    transcripts_df = gdf[gdf.Feature == "transcript"][["Chromosome", "Start", "End", "Strand"]].copy()

    if transcripts_df.empty and ("transcript_id" in gdf.columns):
        spans = (
            gdf[gdf.Feature == "exon"]
            .groupby(["Chromosome", "Strand", "transcript_id"], as_index=False)
            .agg(Start=("Start", "min"), End=("End", "max"))
        )
        transcripts_df = spans[["Chromosome", "Start", "End", "Strand"]].copy()

    bounds = rmats_all_bounds_df[["Chromosome", "Start", "End", "Strand"]].copy()
    bounds["psi"] = np.nan

    seg01 = _disjoint_labels_with_priority(
        transcripts_df=transcripts_df,
        exons_df=exons_df,
        psi_df=bounds,
    )

    seg01 = seg01.rename(columns={"p_exon": "base_exon"})
    return seg01[["Chromosome", "Start", "End", "Strand", "base_exon"]]


def _map_psi_to_output_paths(cell_type: str, manifest_sub: pd.DataFrame, psi_cols, per_sample_out_dir: str):
    mapping = {}
    sub = manifest_sub.copy()

    rid = None
    if "replicate_id" in sub.columns:
        rid = pd.to_numeric(sub["replicate_id"], errors="coerce")

    used_idx = set()
    for psi_col in psi_cols:
        psi_idx = int(psi_col.split("_")[1])
        chosen_idx = None

        if rid is not None:
            matches = sub.index[rid == psi_idx].tolist()
            if matches:
                chosen_idx = matches[0]

        if chosen_idx is None and "bed_path" in sub.columns:
            contains = sub.index[sub["bed_path"].fillna("").str.contains(psi_col, regex=False)].tolist()
            contains = [i for i in contains if i not in used_idx]
            if contains:
                chosen_idx = contains[0]

        if chosen_idx is not None:
            out_path = str(sub.loc[chosen_idx, "bed_path"])
            used_idx.add(chosen_idx)
        else:
            out_path = str(Path(per_sample_out_dir) / f"{cell_type}_soft_{psi_col}.bed")

        mapping[psi_col] = out_path

    return mapping


def _project_values_to_segments(segments_df: pd.DataFrame, value_df: pd.DataFrame, value_col: str, default_value: float = 0.0):
    S = segments_df[["Chromosome", "Start", "End", "Strand"]].copy()
    V = value_df[["Chromosome", "Start", "End", "Strand", value_col]].copy() if not value_df.empty else pd.DataFrame(columns=["Chromosome", "Start", "End", "Strand", value_col])

    S = S.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    V = V.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True) if not V.empty else V

    out = []
    keys = S[["Chromosome", "Strand"]].drop_duplicates().itertuples(index=False, name=None)
    for chrom, strand in keys:
        seg = S[(S.Chromosome == chrom) & (S.Strand == strand)]
        val = V[(V.Chromosome == chrom) & (V.Strand == strand)] if not V.empty else V

        seg_starts = seg["Start"].to_numpy()
        seg_ends = seg["End"].to_numpy()

        val_starts = val["Start"].to_numpy() if val is not None and not val.empty else np.array([], dtype=int)
        val_ends = val["End"].to_numpy() if val is not None and not val.empty else np.array([], dtype=int)
        val_vals = val[value_col].to_numpy() if val is not None and not val.empty else np.array([], dtype=float)

        j = 0
        for i in range(len(seg)):
            s0 = int(seg_starts[i])
            s1 = int(seg_ends[i])

            while j < len(val_ends) and int(val_ends[j]) <= s0:
                j += 1

            k = j
            vals = []
            while k < len(val_starts) and int(val_starts[k]) < s1:
                if int(val_ends[k]) > s0:
                    vals.append(float(val_vals[k]))
                k += 1

            assigned = float(np.mean(vals)) if vals else float(default_value)
            out.append((chrom, s0, s1, strand, assigned))

    return pd.DataFrame(out, columns=["Chromosome", "Start", "End", "Strand", value_col])


def regenerate_per_sample_beds_from_rmats(manifest: pd.DataFrame, args) -> pd.DataFrame:
    if manifest.empty:
        return manifest

    updated = manifest.copy()
    if "bed_path" not in updated.columns:
        updated["bed_path"] = ""

    for cell_type, sub in manifest.groupby("cell_type", sort=False):
        produced_rows, _ = regenerate_one_cell_type(cell_type, args, sub)
        for row in produced_rows:
            sample_idx = row["replicate_id"]
            out_path = row["bed_path"]
            row_mask = (
                (updated["cell_type"] == cell_type)
                & (updated.get("replicate_id", pd.Series([None] * len(updated))).astype(str) == str(sample_idx))
            )
            if row_mask.any():
                updated.loc[row_mask, "bed_path"] = str(out_path)
            else:
                contains_mask = (
                    (updated["cell_type"] == cell_type)
                    & updated["bed_path"].fillna("").str.contains(f"psi_{sample_idx}", regex=False)
                )
                if contains_mask.any():
                    updated.loc[contains_mask, "bed_path"] = str(out_path)

    return updated


def build_manifest_from_cell_types(args) -> pd.DataFrame:
    rows = []
    for cell_type in _parse_cell_types(args.cell_types):
        produced_rows, _ = regenerate_one_cell_type(cell_type, args, sub=None)
        rows.extend(produced_rows)

    cols = ["sample_id", "cell_type", "bed_path", "replicate_id", "rmats_se", "rmats_ri", "rmats_a3ss", "rmats_a5ss", "gtf"]
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows, columns=cols)


def run_master_mode(args):
    """Run in master BED creation mode."""
    if args.regenerate_beds:
        print("[INFO] Master mode will REGENERATE per-sample BEDs from rMATS first, then aggregate masters.")
    else:
        print("[INFO] Master mode aggregates existing BED files listed in manifest (--no-regenerate-beds).")
    manifest = None
    if args.manifest:
        manifest_path = Path(args.manifest)
        if not manifest_path.exists():
            print(f"[ERROR] Manifest file not found: {manifest_path}")
            return 1
        print(f"[INFO] Loading manifest from {manifest_path}")
        try:
            manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
        except Exception as e:
            print(f"[ERROR] reading manifest: {e}")
            return 1
        required_cols = {"cell_type", "bed_path"}
        if not required_cols.issubset(manifest.columns):
            print(f"[ERROR] Manifest must contain columns: {required_cols}")
            print(f"[ERROR] Found columns: {set(manifest.columns)}")
            return 1
    
    # Regenerate per-sample BEDs from rMATS first (default)
    if args.regenerate_beds:
        if manifest is not None:
            manifest = regenerate_per_sample_beds_from_rmats(manifest, args)
        else:
            print("[INFO] No input manifest provided; building manifest from regenerated per-sample BEDs.")
            manifest = build_manifest_from_cell_types(args)

        if manifest.empty:
            print("[ERROR] No regenerated per-sample BEDs were produced. Cannot build master BEDs.")
            return 1

        manifest_out = args.manifest_out
        if not manifest_out:
            manifest_out = str(Path(args.out_dir) / "manifest.generated.tsv")
        out_manifest = Path(manifest_out)
        out_manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.to_csv(out_manifest, sep="\t", index=False)
        print(f"[INFO] Wrote manifest to {out_manifest}")
    else:
        if manifest is None:
            print("[ERROR] --no-regenerate-beds requires --manifest (no BED source provided).")
            return 1

    # Group BED files by cell type
    bed_files = {}
    for ct, sub in manifest.groupby("cell_type", sort=False):
        bed_files[ct] = sub["bed_path"].tolist()
    
    print(f"[INFO] Found {len(bed_files)} cell type(s): {', '.join(sorted(bed_files.keys()))}")
    
    # Validate that BED files exist
    missing_files = []
    for ct, paths in bed_files.items():
        for path in paths:
            if not Path(path).exists():
                missing_files.append(path)
    
    if missing_files:
        print(f"[WARN] {len(missing_files)} BED file(s) not found:")
        for f in missing_files[:5]:
            print(f"[WARN]   {f}")
        if len(missing_files) > 5:
            print(f"[WARN]   ... and {len(missing_files) - 5} more")
    
    # Create master BED files
    create_master_bed_files(bed_files, args.out_dir)
    print("[INFO] Done!")
    return 0


def run_single_mode(args):
    """Run in single BED generation mode."""
    print(f"[INFO] Loading SE from {args.rmats_se}")
    rmats_se = load_rmats_se(args.rmats_se, args.min_coverage)
    
    print(f"[INFO] Loading RI from {args.rmats_ri}")
    rmats_ri = load_rmats_ri(args.rmats_ri, args.min_coverage)

    print(f"[INFO] Loading A3SS from {args.rmats_a3ss}")
    rmats_a3ss = load_rmats_alt_site(args.rmats_a3ss, args.min_coverage)

    print(f"[INFO] Loading A5SS from {args.rmats_a5ss}")
    rmats_a5ss = load_rmats_alt_site(args.rmats_a5ss, args.min_coverage)
    
    print(f"[INFO] Loading GTF from {args.gtf}")
    stringtie_gr = load_stringtie_gtf(args.gtf)

    print("[INFO] Generating Labels...")
    df = generate_soft_labels(
        stringtie_gr,
        rmats_se,
        rmats_ri=rmats_ri,
        rmats_a3ss=rmats_a3ss,
        rmats_a5ss=rmats_a5ss,
    )

    out_path = Path(args.out_bed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Add default confidence if not present
    if 'confidence' not in df.columns:
        df['confidence'] = 0.0
    
    df[['Chromosome', 'Start', 'End', 'p_exon', 'confidence']].to_csv(
        out_path, sep="\t", header=False, index=False
    )

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        df['p_exon'].describe().rename_axis('stat').to_frame('value').to_csv(
            summary_path, sep="\t"
        )

    print(f"[INFO] Wrote {len(df)} regions to {out_path}")
    return 0


def main():
    args = parse_args()
    
    # Determine mode
    if args.manifest or args.out_dir:
        # Master BED creation mode (manifest optional if regenerating)
        if not args.out_dir:
            print("[ERROR] --out-dir is required for master BED mode")
            return 1
        return run_master_mode(args)
    else:
        # Single BED generation mode
        return run_single_mode(args)

if __name__ == "__main__":
    exit(main())