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
    def __init__(self, fasta_path: str):
        if pysam is None:
            raise ImportError("pysam is required. pip install pysam")
        self.fasta = pysam.FastaFile(fasta_path)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        return self.fasta.fetch(chrom, start, end)

def load_label_intervals(bed_path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(
            bed_path,
            sep="\t",
            names=["Chromosome", "Start", "End", "p_exon"],
            dtype={"Chromosome": str},
        )
    except Exception as e:
        print(f"Error reading {bed_path}: {e}")
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
    # Fix: Ensure robust loading by handling NaNs before int conversion
    df = pd.read_csv(path, sep="\t", dtype={"Chromosome": str})
    
    required = {"Chromosome", "Start", "End"}
    if not required.issubset(df.columns):
         # Try to be helpful if columns are missing (e.g. maybe no header?)
         # But usually strictly requiring headers is safer for pipelines.
         pass

    # 1. Coerce bad values to NaN
    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    
    # 2. Drop rows with NaNs
    original_len = len(df)
    df = df.dropna(subset=["Start", "End"]).reset_index(drop=True)
    if len(df) < original_len:
        print(f"Warning: Dropped {original_len - len(df)} transcript rows due to invalid coordinates.")

    # 3. Safe to cast to int
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    
    # 4. Filter logic
    df = df[df["Start"] < df["End"]].reset_index(drop=True)
    
    return df

# ----------------------------
# Optimized Windowing + Rasterization
# ----------------------------

def iter_windows_from_regions(regions_df: pd.DataFrame, window_bp: int, stride_bp: int):
    """Yield genomic windows from transcript regions."""
    for chrom, sub in regions_df.groupby("Chromosome"):
        start = int(sub["Start"].min())
        end = int(sub["End"].max())
        if end - start < window_bp: continue
        for w_start in range(start, end - window_bp + 1, stride_bp):
            yield chrom, int(w_start), int(w_start + window_bp)

def rasterize_window_vectorized(
    window_start: int,
    window_end: int,
    bin_size: int,
    starts_arr: np.ndarray,
    ends_arr: np.ndarray,
    psis_arr: np.ndarray,
    n_bins: int,
):
    """
    Highly optimized version using NumPy arrays instead of Pandas iteration.
    """
    targets = np.zeros(n_bins, dtype=np.float32)
    events = np.zeros(n_bins, dtype=np.float32)

    if len(starts_arr) == 0:
        return targets, events

    # 1. Vectorized Filter: Find only exons overlapping this window
    #    (Start < WindowEnd) AND (End > WindowStart)
    overlap_mask = (starts_arr < window_end) & (ends_arr > window_start)
    
    if not np.any(overlap_mask):
        return targets, events

    # Subset arrays to relevant exons
    s_sub = starts_arr[overlap_mask]
    e_sub = ends_arr[overlap_mask]
    p_sub = psis_arr[overlap_mask]

    # 2. Iterate only over relevant exons (much faster than iterating all)
    #    We stick to a loop here because n_exons_in_window is usually small (10-50),
    #    whereas n_bins is large (16k). Fully vectorizing bins usually consumes too much RAM.
    for s, e, p in zip(s_sub, e_sub, p_sub):
        if np.isnan(p): continue
        
        # Clip to window
        s_in = max(s, window_start)
        e_in = min(e, window_end)
        
        # Calculate bin indices
        b0 = (s_in - window_start) // bin_size
        b1 = (e_in - window_start + bin_size - 1) // bin_size
        
        # Python loop over the bins covered by this specific exon
        # (Usually an exon covers only 3-5 bins, so this is fast)
        for b in range(int(b0), int(b1)):
            bin_s = window_start + b * bin_size
            bin_e = bin_s + bin_size
            
            # Intersection
            ov = max(0, min(e_in, bin_e) - max(s_in, bin_s))
            if ov <= 0: continue

            targets[b] += ov * float(p)
            events[b] += ov

    # Normalize
    nz = events > 0
    targets[nz] /= events[nz]

    return targets, events

def one_hot_encode(seq: str) -> np.ndarray:
    mapping = {
        "A": (1, 0, 0, 0), "C": (0, 1, 0, 0), "G": (0, 0, 1, 0), "T": (0, 0, 0, 1),
        "N": (0, 0, 0, 0), "a": (1, 0, 0, 0), "c": (0, 1, 0, 0), "g": (0, 0, 1, 0), "t": (0, 0, 0, 1)
    }
    out = np.zeros((len(seq), 4), dtype=np.float32)
    for i, b in enumerate(seq):
        out[i] = mapping.get(b, (0, 0, 0, 0))
    return out

# ----------------------------
# Sharding & Writing
# ----------------------------

def save_shard(examples: List[dict], out_dir: str, split: str, shard_id: int):
    out = Path(out_dir) / split
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"shard_{shard_id:04d}.npz"
    np.savez_compressed(
        path,
        sequence=np.stack([e["sequence"] for e in examples]),
        targets=np.stack([e["targets"] for e in examples]),
        mask=np.stack([e["mask"] for e in examples]),
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
    shard_prefix: int = 0
):
    if window_bp % bin_size != 0: raise ValueError("Window/Bin mismatch")
    n_bins = window_bp // bin_size
    n_tasks = len(celltype_to_task_idx)
    
    # 1. Pre-load and NumPy-ify Data
    print(f"Pre-loading BED files for {n_tasks} cell types...")
    # data_cache[ct][chrom] = (starts_arr, ends_arr, psis_arr)
    data_cache = {ct: {} for ct in bed_files_by_celltype}
    
    for ct, bed_paths in bed_files_by_celltype.items():
        for bpath in bed_paths:
            df = load_label_intervals(bpath)
            if df.empty: continue
            for chrom, group in df.groupby("Chromosome"):
                if chrom not in data_cache[ct]:
                    data_cache[ct][chrom] = []
                # Store as tuple of numpy arrays for speed
                data_cache[ct][chrom].append((
                    group["Start"].values, 
                    group["End"].values, 
                    group["p_exon"].values
                ))

    # 2. Window Generation
    FLANK = 262_144
    regions_df = transcript_df.copy()
    regions_df["Start"] = (regions_df["Start"] - FLANK).clip(lower=0)
    regions_df["End"] = regions_df["End"] + FLANK
    regions_df = regions_df[["Chromosome","Start","End"]]

    windows = list(iter_windows_from_regions(regions_df, window_bp, stride_bp))
    print(f"Total windows to scan: {len(windows)}")
    
    examples: List[dict] = []
    shard_id = shard_prefix
    emitted = 0
    
    n_signal = 0
    n_noise = 0

    for chrom, w_start, w_end in windows:
        seq = fasta.fetch(chrom, w_start, w_end)
        if len(seq) != window_bp: continue
        if seq.upper().count("N") / window_bp > 0.5: continue
        
        y = np.zeros((n_bins, n_tasks), dtype=np.float32)
        mask = np.zeros((n_bins, n_tasks), dtype=np.float32)
        window_has_signal = False
        
        for ct, task_idx in celltype_to_task_idx.items():
            replicates_data = data_cache[ct].get(chrom, [])
            
            accum_targets = np.zeros(n_bins, dtype=np.float32)
            accum_mask = np.zeros(n_bins, dtype=np.float32)
            
            for (s_arr, e_arr, p_arr) in replicates_data:
                # Optimized rasterization call
                targ, events = rasterize_window_vectorized(
                    w_start, w_end, bin_size, s_arr, e_arr, p_arr, n_bins
                )
                
                if np.any(events > 0):
                    rep_mask = (events > 0).astype(np.float32)
                    accum_targets += targ 
                    accum_mask += rep_mask

            # Average across replicates
            with np.errstate(divide='ignore', invalid='ignore'):
                y[:, task_idx] = np.divide(accum_targets, accum_mask)
                y[accum_mask == 0, task_idx] = 0.0
            
            mask[:, task_idx] = (accum_mask > 0).astype(np.float32)
        
        if np.any(mask > 0):
            window_has_signal = True

        # Balancing Logic
        include = True
        if subset_mode == "event_only":
            include = window_has_signal
        elif subset_mode == "balanced":
            if window_has_signal:
                n_signal += 1
            else:
                if n_noise < n_signal:
                    n_noise += 1
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
    
    print(f"Done. Emitted {emitted} examples.")

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--transcripts", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--split", default="train")
    ap.add_argument("--window-bp", type=int, default=524288)
    ap.add_argument("--bin-size", type=int, default=32)
    ap.add_argument("--stride-bp", type=int, default=262144)
    ap.add_argument("--shard-size", type=int, default=256)
    ap.add_argument("--max-examples", type=int, default=None)
    ap.add_argument("--subset-mode", choices=["all", "event_only", "balanced"], default="all")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--chroms", type=str, default=None)
    ap.add_argument("--shard-prefix", type=int, default=0, help="Start numbering shards from here (for parallelism)")
    return ap.parse_args()

def main():
    args = parse_args()
    fasta = FastaReader(args.fasta)
    transcript_df = load_transcripts_tsv(args.transcripts)
    
    if args.chroms:
        keep = set(c.strip() for c in args.chroms.split(","))
        transcript_df = transcript_df[transcript_df["Chromosome"].isin(keep)]

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str)
    bed_files = {}
    for ct, sub in manifest.groupby("cell_type", sort=False):
        bed_files[ct] = sub["bed_path"].tolist()
    
    celltypes = sorted(bed_files.keys())
    task_map = {ct: i for i, ct in enumerate(celltypes)}

    write_dataset_shards(
        fasta=fasta, transcript_df=transcript_df, bed_files_by_celltype=bed_files,
        celltype_to_task_idx=task_map, out_dir=args.out_dir, split=args.split,
        window_bp=args.window_bp, bin_size=args.bin_size, stride_bp=args.stride_bp,
        shard_size=args.shard_size, max_examples=args.max_examples,
        subset_mode=args.subset_mode, seed=args.seed, shard_prefix=args.shard_prefix
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