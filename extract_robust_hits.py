#!/usr/bin/env python3
"""
Extract 'Robust' Splicing Events from rMATS results.

This script filters rMATS output to find events that are both:
1. Statistically Significant (low FDR)
2. Biologically Relevant (high magnitude of change / dPSI)

Usage:
    python extract_robust_events.py <rMATS_result_dir> [options]
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Default thresholds
DEFAULT_FDR = 0.05       # Significance cutoff
DEFAULT_DPSI = 0.1       # 10% splicing difference cutoff (Effect Size)
EVENT_TYPES = ["SE", "RI", "MXE", "A3SS", "A5SS"]

def parse_args():
    parser = argparse.ArgumentParser(description="Filter rMATS results for robust biological events.")
    parser.add_argument("result_dir", type=Path, help="Directory containing rMATS text files")
    parser.add_argument("--fdr", type=float, default=DEFAULT_FDR, help=f"Max FDR (default: {DEFAULT_FDR})")
    parser.add_argument("--dpsi", type=float, default=DEFAULT_DPSI, help=f"Min abs(IncLevelDifference) (default: {DEFAULT_DPSI})")
    parser.add_argument("--top", type=int, default=50, help="Number of top events to display in summary")
    parser.add_argument("--output", type=str, default="robust_events_summary.csv", help="Output filename for combined hits")
    return parser.parse_args()

def process_file(filepath, fdr_cutoff, dpsi_cutoff):
    if not filepath.exists():
        return None
    
    try:
        df = pd.read_csv(filepath, sep="\t")
    except Exception as e:
        print(f"[WARN] Could not read {filepath.name}: {e}")
        return None

    if df.empty:
        return None

    # rMATS columns of interest
    required_cols = ["GeneID", "geneSymbol", "FDR", "IncLevelDifference", "PValue"]
    
    # Check if columns exist
    if not all(col in df.columns for col in required_cols):
        return None

    # Clean data
    df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")
    df["IncLevelDifference"] = pd.to_numeric(df["IncLevelDifference"], errors="coerce")
    df = df.dropna(subset=["FDR", "IncLevelDifference"])

    # Filter
    # 1. Significance Filter
    mask_sig = df["FDR"] <= fdr_cutoff
    
    # 2. Biological Magnitude Filter (Absolute Delta PSI)
    mask_bio = df["IncLevelDifference"].abs() >= dpsi_cutoff

    robust_df = df[mask_sig & mask_bio].copy()
    
    # Add helper column for sorting
    robust_df["abs_dPSI"] = robust_df["IncLevelDifference"].abs()
    
    # Sort by magnitude of change (biggest switches first)
    robust_df = robust_df.sort_values(by="abs_dPSI", ascending=False)
    
    return robust_df

def main():
    args = parse_args()
    
    if not args.result_dir.exists():
        print(f"Error: Directory '{args.result_dir}' not found.")
        sys.exit(1)

    print(f"Scanning {args.result_dir}...")
    print(f"Criteria: FDR <= {args.fdr}, |dPSI| >= {args.dpsi}")
    
    all_hits = []

    for event in EVENT_TYPES:
        filename = f"{event}.MATS.JC.txt"
        filepath = args.result_dir / filename
        
        result = process_file(filepath, args.fdr, args.dpsi)
        
        if result is not None and not result.empty:
            count = len(result)
            print(f"  -> {event}: Found {count} robust events.")
            
            # Add event type column for the combined report
            result["EventType"] = event
            all_hits.append(result)
        else:
            print(f"  -> {event}: No robust events found.")

    if all_hits:
        combined = pd.concat(all_hits, ignore_index=True)
        # Sort combined list by absolute dPSI
        combined = combined.sort_values(by="abs_dPSI", ascending=False)
        
        # Save to CSV
        out_path = args.result_dir / args.output
        combined.to_csv(out_path, index=False)
        print(f"\n[OK] Saved {len(combined)} robust events to: {out_path}")
        
        # Preview top hits
        print("\n=== Top 10 Most Robust Switching Events ===")
        print(f"{'Gene':<15} {'Type':<6} {'dPSI':<8} {'FDR':<10}")
        print("-" * 45)
        for _, row in combined.head(10).iterrows():
            gene = str(row['geneSymbol'])[:14]
            etype = row['EventType']
            dpsi = f"{row['IncLevelDifference']:.3f}"
            fdr = f"{row['FDR']:.2e}"
            print(f"{gene:<15} {etype:<6} {dpsi:<8} {fdr:<10}")
            
    else:
        print("\n[!] No events met the criteria. Try lowering the dPSI threshold.")

if __name__ == "__main__":
    main()