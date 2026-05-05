#!/usr/bin/env python3
"""
Quickly scan BAM headers to estimate read length and check for compatibility issues
before running rMATS.
Requires: pysam (pip install pysam)
"""

import sys
import pysam
import statistics
from pathlib import Path

def check_bam(bam_path, num_reads=10000):
    path = Path(bam_path)
    if not path.exists():
        print(f"[FAIL] File not found: {bam_path}")
        return

    try:
        # Check index existence (often causes rMATS issues if missing)
        if not (path.with_suffix(".bam.bai").exists() or path.with_suffix(".bai").exists()):
             print(f"[WARN] Index missing for {path.name} (rMATS might complain)")

        with pysam.AlignmentFile(path, "rb") as bam:
            lengths = []
            
            # fast inspection of first N reads
            for i, read in enumerate(bam.fetch()):
                if i >= num_reads:
                    break
                lengths.append(read.query_length)

            if not lengths:
                print(f"[FAIL] {path.name}: No reads found.")
                return

            avg_len = statistics.mean(lengths)
            mode_len = statistics.mode(lengths)
            
            print(f"[OK] {path.name}")
            print(f"     - Read Length (Mode): {mode_len} bp")
            print(f"     - Read Length (Avg):  {avg_len:.1f} bp")
            
            # Simple check for strandedness (heuristic)
            # This is hard to do perfectly without a reference, but we can check flags
            # to see if it's paired properly.
            
    except Exception as e:
        print(f"[ERR] {path.name}: {e}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python check_rmats_headers.py <list_of_bams.txt> [or single .bam]")
        sys.exit(1)

    targets = []
    # If arg is a text file list (like B_Cell_bams.txt)
    arg_path = Path(sys.argv[1])
    
    if arg_path.suffix == ".txt":
        with open(arg_path) as f:
            # Handle comma separated lists or newlines
            content = f.read().strip()
            if "," in content:
                targets = content.split(",")
            else:
                targets = content.splitlines()
    else:
        targets = sys.argv[1:]

    print(f"Checking {len(targets)} BAM files...")
    print("-" * 40)
    for t in targets:
        check_bam(t.strip())

if __name__ == "__main__":
    main()