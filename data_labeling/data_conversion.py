#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import time
import logging

import numpy as np
import pandas as pd

try:
    import pysam
except ImportError as e:
    pysam = None

# set logging level
logging.basicConfig(level=logging.INFO)

# ----------------------------
# IO helpers
# ----------------------------

class FastaReader:
    """Thin wrapper around pysam.FastaFile with 0-based half-open fetch()."""
    def __init__(self, fasta_path: str):
        if pysam is None:
            raise ImportError("pysam is required for FASTA access. Install with: pip install pysam")
        self.fasta = pysam.FastaFile(fasta_path)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        return self.fasta.fetch(chrom, start, end)

def load_label_intervals(bed_path: str) -> pd.DataFrame:
    """BED must be 0-based, half-open, with 4 columns: chrom start end p_exon."""
    try:
        df = pd.read_csv(
            bed_path,
            sep="\t",
            names=["Chromosome", "Start", "End", "p_exon"],
            dtype={"Chromosome": str},
        )
    except Exception as e:
        logging.error(f"Error reading {bed_path}: {e}")
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "p_exon"])

    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df = df.dropna(subset=["Start", "End"]).reset_index(drop=True)
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["p_exon"] = df["p_exon"].astype(float)
    df = df[df["Start"] < df["End"]].reset_index(drop=True)
    return df

def load_transcripts_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"Chromosome": str})
    required = {"Chromosome", "Start", "End"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Transcript TSV missing columns: {sorted(missing)}")

    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df = df.dropna(subset=["Start", "End"]).reset_index(drop=True)
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return df

# ----------------------------
# Windowing + Rasterization
# ----------------------------

def iter_windows_from_regions(
    regions_df: pd.DataFrame,
    window_bp: int,
    stride_bp: int,
):
    """Yield genomic windows from the UNION of transcript regions per chromosome."""
    required = {"Chromosome", "Start", "End"}
    missing = required - set(regions_df.columns)
    if missing:
        raise ValueError(f"regions_df missing columns: {sorted(missing)}")

    for chrom, sub in regions_df.groupby("Chromosome"):
        start = int(sub["Start"].min())
        end = int(sub["End"].max())

        if end - start < window_bp:
            continue

        for w_start in range(start, end - window_bp + 1, stride_bp):
            yield chrom, int(w_start), int(w_start + window_bp)

def rasterize_window(
    window_start: int,
    window_end: int,
    bin_size: int,
    label_df: pd.DataFrame,
    n_bins: int,
):
    """
    Returns:
      targets[b] = mean p_exon (weighted by overlap)
      events[b]  = amount of overlap with valid events
    """
    targets = np.zeros(n_bins, dtype=np.float32)
    events = np.zeros(n_bins, dtype=np.float32)

    if label_df.empty:
        return targets, events

    # Optimization: Filter roughly before iterating
    # (Though logic below handles boundaries correctly, this speeds up large DFs)
    # subset = label_df[(label_df["End"] > window_start) & (label_df["Start"] < window_end)]
    
    # Actually, because we pre-filter by chromosome in the main loop, 
    # and label_df passed here is already specific, we can just iterate.
    # For massive BEDs, using searchsorted on Start/End arrays is faster, 
    # but basic iteration is okay for PSI events which are usually sparse.
    
    for _, row in label_df.iterrows():
        s = int(row.Start)
        e = int(row.End)
        
        # Skip if no overlap
        if e <= window_start or s >= window_end:
            continue

        # Constrain to window
        s_in = max(s, window_start)
        e_in = min(e, window_end)
        
        b0 = (s_in - window_start) // bin_size
        b1 = (e_in - window_start + bin_size - 1) // bin_size

        p = row.p_exon
        has_event = not np.isnan(p)

        if not has_event:
            continue

        for b in range(int(b0), int(b1)):
            bin_s = window_start + b * bin_size
            bin_e = bin_s + bin_size
            
            # Intersection of bin and exon segment
            ov = max(0, min(e_in, bin_e) - max(s_in, bin_s))
            if ov <= 0:
                continue

            targets[b] += ov * float(p)
            events[b] += ov

    # Normalize by coverage amount
    nz = events > 0
    targets[nz] /= events[nz]

    return targets, events

def one_hot_encode(seq: str) -> np.ndarray:
    mapping = {
        "A": (1, 0, 0, 0), "C": (0, 1, 0, 0), "G": (0, 0, 1, 0), "T": (0, 0, 0, 1),
        "N": (0, 0, 0, 0), "a": (1, 0, 0, 0), "c": (0, 1, 0, 0), "g": (0, 0, 1, 0), "t": (0, 0, 0, 1)
    }
    # Pre-allocate
    out = np.zeros((len(seq), 4), dtype=np.float32)
    for i, b in enumerate(seq):
        out[i] = mapping.get(b, (0, 0, 0, 0))
    return out

# ----------------------------
# Dataset Generation
# ----------------------------

def save_shard(examples: List[dict], out_dir: str, split: str, shard_id: int):
    out = Path(out_dir) / split
    out.mkdir(parents=True, exist_ok=True)
    
    path = out / f"shard_{shard_id:04d}.npz"
    
    # Stack lists into arrays
    np.savez_compressed(
        path,
        sequence=np.stack([e["sequence"] for e in examples]),
        targets=np.stack([e["targets"] for e in examples]),
        mask=np.stack([e["mask"] for e in examples]),
        # Note: cell_type metadata is less relevant per-example in multitask 
        # because every example contains ALL cell types. 
        # We store the list of cell types in a separate metadata file usually.
        chrom=np.array([e["chrom"] for e in examples], dtype="U"),
        window_start=np.array([e["window_start"] for e in examples], dtype=np.int64),
        window_end=np.array([e["window_end"] for e in examples], dtype=np.int64),
        has_signal=np.array([e["has_signal"] for e in examples], dtype=np.float32),
    )

def write_dataset_shards(
    fasta: FastaReader,
    transcript_df: pd.DataFrame,
    bed_files_by_celltype: Dict[str, List[str]],
    celltype_to_task_idx: Dict[str, int],
    out_dir: str,
    split: str,
    window_bp: int = 524288,
    bin_size: int = 32,
    stride_bp: int = 262144,
    shard_size: int = 128,
    max_examples: Optional[int] = None,
    subset_mode: str = "all",
    seed: int = 42,
):
    """
    Revised logic:
    1. Pre-load ALL bed files into memory (grouped by CellType -> Chromosome).
    2. Iterate Windows.
    3. For each window, rasterize ALL cell types (averaging replicates).
    4. Save Dense Multitask Example.
    """
    
    if window_bp % bin_size != 0:
        raise ValueError(f"window_bp ({window_bp}) must be divisible by bin_size ({bin_size})")

    n_bins = window_bp // bin_size
    n_tasks = len(celltype_to_task_idx)
    
    # --- 1. Pre-load Data (Crucial for Speed) ---
    logging.info(f"Pre-loading BED files for {n_tasks} cell types...")
    # Structure: data_cache[cell_type][chrom] = list of DataFrames (one per replicate)
    data_cache = {ct: {} for ct in bed_files_by_celltype}
    
    for ct, bed_paths in bed_files_by_celltype.items():
        logging.info(f"  Loading {len(bed_paths)} replicates for {ct}...")
        for bpath in bed_paths:
            df = load_label_intervals(bpath)
            if df.empty: continue
            
            # Group by chrom for fast access
            for chrom, group in df.groupby("Chromosome"):
                if chrom not in data_cache[ct]:
                    data_cache[ct][chrom] = []
                data_cache[ct][chrom].append(group)

    # --- 2. Window Generation ---
    FLANK = 262_144
    regions_df = transcript_df.copy()
    regions_df["Start"] = (regions_df["Start"] - FLANK).clip(lower=0)
    regions_df["End"] = regions_df["End"] + FLANK
    regions_df = regions_df[["Chromosome","Start","End"]]

    windows = list(iter_windows_from_regions(regions_df, window_bp, stride_bp))
    logging.info(f"Total windows to scan: {len(windows)}")
    
    examples: List[dict] = []
    shard_id = 0
    emitted = 0
    
    rng = np.random.default_rng(seed)
    # Counters for balancing (Tracking "Any Signal" vs "No Signal")
    n_signal = 0
    n_noise = 0

    # --- 3. Main Loop ---
    for chrom, w_start, w_end in windows:
        seq = fasta.fetch(chrom, w_start, w_end)
        if len(seq) != window_bp: continue
        if seq.upper().count("N") / window_bp > 0.5: continue
        
        # Prepare Dense Target Matrix
        # Shape: (Bins, Tasks)
        y = np.zeros((n_bins, n_tasks), dtype=np.float32)
        mask = np.zeros((n_bins, n_tasks), dtype=np.float32)
        
        window_has_signal = False
        
        # Fill each task column
        for ct, task_idx in celltype_to_task_idx.items():
            replicates_in_chrom = data_cache[ct].get(chrom, [])
            
            # We will average the rasterized targets across replicates
            accum_targets = np.zeros(n_bins, dtype=np.float32)
            accum_mask = np.zeros(n_bins, dtype=np.float32)
            valid_reps = 0
            
            for rep_df in replicates_in_chrom:
                # Fast filter: Check if any events in this DF overlap window
                # (Simple min/max check on the whole DF group is fast)
                # But since we already grouped by chrom, we just rasterize.
                # rasterize_window handles the binary search/overlap logic.
                
                targ, events = rasterize_window(w_start, w_end, bin_size, rep_df, n_bins)
                
                # If this replicate had any events in this window
                if np.any(events > 0):
                    rep_mask = (events > 0).astype(np.float32)
                    accum_targets += targ  # targ is already mean PSI for that replicate
                    accum_mask += rep_mask
                    valid_reps += 1

            if valid_reps > 0:
                # Average the PSI values across the replicates that had coverage
                # Note: We divide by valid_reps? Or should we average only where mask overlaps?
                # Simpler robust approach:
                # y = Sum(Targets) / Valid_Reps. 
                # This assumes replicates cover similar regions.
                
                # More precise: Average only non-zero contributions? 
                # No, 0 is a valid PSI (intron).
                # We need to know if the replicate *covered* the region.
                # Our rasterize returns "events" which is coverage.
                
                # Current logic: If a bin was covered by Replicate A but not Replicate B,
                # Replicate B contributes nothing.
                # We should average by the number of replicates that covered that specific bin.
                
                # Let's refine the accumulation:
                # We need per-bin counts of how many replicates covered it.
                # accum_mask effectively holds this count!
                
                with np.errstate(divide='ignore', invalid='ignore'):
                    y[:, task_idx] = np.divide(accum_targets, accum_mask)
                    # Where accum_mask is 0, result is NaN/Inf. Set to 0.
                    y[accum_mask == 0, task_idx] = 0.0
                
                # Final mask is 1 if ANY replicate had coverage
                mask[:, task_idx] = (accum_mask > 0).astype(np.float32)
        
        # Check overall signal
        if np.any(mask > 0):
            window_has_signal = True

        # --- 4. Filtering / Balancing ---
        include = True
        
        if subset_mode == "event_only":
            include = window_has_signal
            
        elif subset_mode == "balanced":
            # Simple global balancing: Keep 1:1 ratio of Signal vs No-Signal windows
            if window_has_signal:
                n_signal += 1
                include = True
            else:
                if n_noise < n_signal:
                    n_noise += 1
                    include = True
                else:
                    include = False
                    
        if include:
            examples.append({
                "sequence": one_hot_encode(seq),
                "targets": y,
                "mask": mask,
                "chrom": chrom,
                "window_start": int(w_start),
                "window_end": int(w_end),
                "has_signal": float(window_has_signal),
            })
            
            emitted += 1
            if max_examples is not None and emitted >= max_examples:
                save_shard(examples, out_dir, split, shard_id)
                return

            if len(examples) >= shard_size:
                save_shard(examples, out_dir, split, shard_id)
                shard_id += 1
                examples = []

    if examples:
        save_shard(examples, out_dir, split, shard_id)
    
    logging.info(f"Done. Emitted {emitted} examples.")


# ----------------------------
# CLI
# ----------------------------

def parse_args():
    ap = argparse.ArgumentParser(description="Convert per-sample BED soft labels into DENSE MULTITASK windowed shards.")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--transcripts", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--split", default="train")
    ap.add_argument("--window-bp", type=int, default=131072)
    ap.add_argument("--bin-size", type=int, default=32)
    ap.add_argument("--stride-bp", type=int, default=65536)
    ap.add_argument("--shard-size", type=int, default=128)
    ap.add_argument("--max-examples", type=int, default=None)
    ap.add_argument("--subset-mode", choices=["all", "event_only", "balanced"], default="all")
    ap.add_argument("--seed", type=int, default=42)
    
    # Manifest is mandatory now
    ap.add_argument("--manifest", required=True, help="TSV mapping cell_type to bed_path.")
    ap.add_argument("--chroms", type=str, default=None)
    
    # Compatibility args (ignored)
    ap.add_argument("--cell-type", action="append", default=[], help="Ignored, use manifest")
    ap.add_argument("--bed", action="append", default=[], help="Ignored, use manifest")
    
    return ap.parse_args()

def main():
    args = parse_args()

    fasta = FastaReader(args.fasta)
    transcript_df = load_transcripts_tsv(args.transcripts)

    if args.chroms is not None:
        keep = set(c.strip() for c in args.chroms.split(","))
        transcript_df = transcript_df[transcript_df["Chromosome"].isin(keep)]

    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found at {manifest_path}")

    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    
    # Group beds by cell type
    bed_files_by_celltype: Dict[str, List[str]] = {}
    for ct, sub in manifest.groupby("cell_type", sort=False):
        bed_files_by_celltype[ct] = sub["bed_path"].tolist()

    # Define task order deterministically
    celltypes = sorted(bed_files_by_celltype.keys())
    celltype_to_task_idx = {ct: i for i, ct in enumerate(celltypes)}

    logging.info("Building DENSE multitask dataset.")
    logging.info("Task Mapping:")
    for ct, idx in celltype_to_task_idx.items():
        logging.info(f"  Column {idx}: {ct} ({len(bed_files_by_celltype[ct])} replicates)")

    write_dataset_shards(
        fasta=fasta,
        transcript_df=transcript_df,
        bed_files_by_celltype=bed_files_by_celltype,
        celltype_to_task_idx=celltype_to_task_idx,
        out_dir=args.out_dir,
        split=args.split,
        window_bp=args.window_bp,
        bin_size=args.bin_size,
        stride_bp=args.stride_bp,
        shard_size=args.shard_size,
        max_examples=args.max_examples,
        subset_mode=args.subset_mode,
        seed=args.seed,
    )

if __name__ == "__main__":
    # prepare to record time
    start_time = time.time()

    # execute
    main()

    # compute and log runtime
    end_time = time.time()
    elapsed = end_time - start_time
    hours, rem = divmod(elapsed, 3600)
    minutes, seconds = divmod(rem, 60)
    logging.info(f"Total execution time: {int(hours)}h {int(minutes)}m {seconds:.2f}s")