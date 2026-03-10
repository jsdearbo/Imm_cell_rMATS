#!/usr/bin/env python3
"""
Standalone script to create master BED files by aggregating replicates per cell type.
Averages PSI and confidence values across replicates for matching genomic coordinates.
"""
import argparse
from pathlib import Path
import logging

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def load_label_intervals(bed_path: str) -> pd.DataFrame:
    """Load a BED file with PSI and confidence values."""
    try:
        df = pd.read_csv(
            bed_path,
            sep="\t",
            names=["Chromosome", "Start", "End", "PSI", "confidence_value"],
            dtype={"Chromosome": str},
        )
    except Exception as e:
        logging.error(f"Error reading {bed_path}: {e}")
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "PSI", "confidence_value"])

    # If confidence_value is all NaN (e.g. file had only 4 cols), fill with 0
    if df["confidence_value"].isna().all():
        df["confidence_value"] = 0.0
    
    # Fill remaining NaNs (if any mixed rows)
    df["confidence_value"] = df["confidence_value"].fillna(0.0)

    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df = df.dropna(subset=["Start", "End"]).reset_index(drop=True)
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["PSI"] = df["PSI"].astype(float)
    df["confidence_value"] = df["confidence_value"].astype(float)
    df = df[df["Start"] < df["End"]].reset_index(drop=True)
    return df


def create_master_bed_files(bed_files_by_celltype: dict, out_dir: str):
    """
    Creates aggregated master BED files for each cell type.
    Averages PSI and confidence_value across replicates for overlapping intervals.
    """
    master_dir = Path(out_dir) / "master_beds"
    master_dir.mkdir(parents=True, exist_ok=True)
    
    logging.info(f"Creating master BED files in {master_dir}")
    
    for ct, bed_paths in bed_files_by_celltype.items():
        logging.info(f"  Processing {ct}: {len(bed_paths)} replicate(s)")
        
        # Load all intervals from all replicates
        all_intervals = []
        for bpath in bed_paths:
            df = load_label_intervals(bpath)
            if df.empty:
                logging.warning(f"    Skipping empty BED file: {bpath}")
                continue
            all_intervals.append(df)
        
        if not all_intervals:
            logging.warning(f"  No valid intervals for {ct}, skipping.")
            continue
        
        # Concatenate all intervals
        combined = pd.concat(all_intervals, ignore_index=True)
        logging.info(f"    Loaded {len(combined)} total intervals from all replicates")
        
        # Sort by chromosome and position
        combined = combined.sort_values(["Chromosome", "Start", "End"]).reset_index(drop=True)
        
        # Group by exact genomic position and average
        grouped = combined.groupby(["Chromosome", "Start", "End"], as_index=False).agg({
            "PSI": "mean",
            "confidence_value": "mean"
        })
        
        # Write master BED file
        output_path = master_dir / f"{ct}_master.bed"
        grouped.to_csv(
            output_path,
            sep="\t",
            header=False,
            index=False,
            columns=["Chromosome", "Start", "End", "PSI", "confidence_value"]
        )
        
        logging.info(f"    Wrote {len(grouped)} unique intervals to {output_path.name}")
    
    logging.info("Master BED file generation complete.")


def parse_args():
    ap = argparse.ArgumentParser(
        description="Create master BED files by aggregating replicates per cell type",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python create_master_beds.py --manifest manifest.tsv --out-dir /path/to/dataset

Manifest file format (tab-separated):
  cell_type    bed_path
  B_cell       /path/to/B_cell_rep1.bed
  B_cell       /path/to/B_cell_rep2.bed
  T_cell       /path/to/T_cell_rep1.bed
        """
    )
    ap.add_argument("--manifest", required=True, help="TSV file with columns: cell_type, bed_path")
    ap.add_argument("--out-dir", required=True, help="Output directory (master_beds/ will be created inside)")
    return ap.parse_args()


def main():
    args = parse_args()
    
    # Validate paths
    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        logging.error(f"Manifest file not found: {args.manifest}")
        return 1
    
    out_dir = Path(args.out_dir)
    
    # Load manifest
    logging.info(f"Loading manifest from {args.manifest}")
    try:
        manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    except Exception as e:
        logging.error(f"Error reading manifest: {e}")
        return 1
    
    # Validate manifest columns
    required_cols = {"cell_type", "bed_path"}
    if not required_cols.issubset(manifest.columns):
        logging.error(f"Manifest must contain columns: {required_cols}")
        logging.error(f"Found columns: {set(manifest.columns)}")
        return 1
    
    # Group BED files by cell type
    bed_files = {}
    for ct, sub in manifest.groupby("cell_type", sort=False):
        bed_files[ct] = sub["bed_path"].tolist()
    
    logging.info(f"Found {len(bed_files)} cell type(s): {', '.join(sorted(bed_files.keys()))}")
    
    # Validate that BED files exist
    missing_files = []
    for ct, paths in bed_files.items():
        for path in paths:
            if not Path(path).exists():
                missing_files.append(path)
    
    if missing_files:
        logging.warning(f"Warning: {len(missing_files)} BED file(s) not found:")
        for f in missing_files[:5]:  # Show first 5
            logging.warning(f"  {f}")
        if len(missing_files) > 5:
            logging.warning(f"  ... and {len(missing_files) - 5} more")
    
    # Create master BED files
    create_master_bed_files(bed_files, args.out_dir)
    
    logging.info("Done!")
    return 0


if __name__ == "__main__":
    exit(main())
