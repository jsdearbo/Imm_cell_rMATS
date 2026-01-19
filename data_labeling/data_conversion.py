#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

try:
    import pysam
except ImportError as e:
    pysam = None


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
        # pysam uses 0-based half-open coordinates already
        return self.fasta.fetch(chrom, start, end)


def load_label_intervals(bed_path: str) -> pd.DataFrame:
    """BED must be 0-based, half-open, with 4 columns: chrom start end p_exon."""
    df = pd.read_csv(
        bed_path,
        sep="\t",
        names=["Chromosome", "Start", "End", "p_exon"],
        dtype={"Chromosome": str},
    )
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["p_exon"] = df["p_exon"].astype(float)
    df = df[df["Start"] < df["End"]].reset_index(drop=True)
    return df


def load_transcripts_tsv(path: str) -> pd.DataFrame:
    """
    Expected columns: Chromosome, Start, End
    Optional: Strand, transcript_id, gene_id, etc.
    Coordinates must be 0-based, half-open.
    """
    df = pd.read_csv(path, sep="\t", dtype={"Chromosome": str})
    required = {"Chromosome", "Start", "End"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Transcript TSV missing columns: {sorted(missing)}")

    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df = df[df["Start"] < df["End"]].reset_index(drop=True)
    return df


# ----------------------------
# Windowing + rasterization
# ----------------------------

# NOTE: Strand is intentionally ignored.
# Borzoi learns strand via RC augmentation.

def iter_windows(transcript_df: pd.DataFrame, window_bp: int, stride_bp: int):
    """
    Yields (chrom, window_start, window_end) for each transcript interval.
    """
    # Assumes transcript_df already filtered for a single chromosome if desired
    for _, row in transcript_df.iterrows():
        start = int(row.Start)
        end = int(row.End)
        chrom = row.Chromosome
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
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Length-weighted mean p_exon per bin. Returns (targets, counts).
    counts[b] is total overlapped bp in that bin (0..bin_size).
    """
    targets = np.zeros(n_bins, dtype=np.float32)
    counts = np.zeros(n_bins, dtype=np.float32)

    # Fast reject if empty
    if label_df.empty:
        return targets, counts

    for _, row in label_df.iterrows():
        s = max(int(row.Start), window_start)
        e = min(int(row.End), window_end)
        if s >= e:
            continue

        # Bin indices (0-based)
        b0 = (s - window_start) // bin_size
        b1 = (e - window_start + bin_size - 1) // bin_size  # ceil

        p = float(row.p_exon)

        for b in range(int(b0), int(b1)):
            bin_s = window_start + b * bin_size
            bin_e = bin_s + bin_size
            ov = max(0, min(e, bin_e) - max(s, bin_s))
            if ov > 0:
                targets[b] += ov * p
                counts[b] += ov

    nz = counts > 0
    targets[nz] /= counts[nz]
    return targets, counts


def one_hot_encode(seq: str) -> np.ndarray:
    mapping = {
        "A": (1, 0, 0, 0),
        "C": (0, 1, 0, 0),
        "G": (0, 0, 1, 0),
        "T": (0, 0, 0, 1),
        "N": (0, 0, 0, 0),
    }
    # Uppercase once
    seq = seq.upper()
    out = np.zeros((len(seq), 4), dtype=np.float32)
    for i, b in enumerate(seq):
        out[i] = mapping.get(b, (0, 0, 0, 0))
    return out


# ----------------------------
# Sharding
# ----------------------------

def save_shard(examples: List[dict], out_dir: str, split: str, shard_id: int):
    out = Path(out_dir) / split
    out.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        out / f"shard_{shard_id:04d}.npz",
        sequence=np.stack([e["sequence"] for e in examples]),
        targets=np.stack([e["targets"] for e in examples]),
        mask=np.stack([e["mask"] for e in examples]),
        cell_type=np.array([e["cell_type"] for e in examples], dtype="U"),
        bed_path=np.array([e["bed_path"] for e in examples], dtype="U"),
        chrom=np.array([e["chrom"] for e in examples], dtype="U"),
        window_start=np.array([e["window_start"] for e in examples], dtype=np.int64),
        window_end=np.array([e["window_end"] for e in examples], dtype=np.int64),
    )


def write_dataset_shards(
    fasta: FastaReader,
    transcript_df: pd.DataFrame,
    bed_files_by_celltype: Dict[str, List[str]],
    celltype_to_task_idx: Dict[str, int],
    out_dir: str,
    split: str,
    window_bp: int = 131072,
    bin_size: int = 32,
    stride_bp: int = 65536,
    shard_size: int = 128,
    max_examples: Optional[int] = None,
):
    if window_bp % bin_size != 0:
        raise ValueError(f"window_bp ({window_bp}) must be divisible by bin_size ({bin_size})")
    n_bins = window_bp // bin_size
    n_tasks = len(celltype_to_task_idx)

    # Basic transcript sanity
    required = {"Chromosome", "Start", "End"}
    missing = required - set(transcript_df.columns)
    if missing:
        raise ValueError(f"Transcript DF missing columns: {sorted(missing)}")

    examples: List[dict] = []
    shard_id = 0
    emitted = 0

    # Iterate cell types (one head per cell type)
    for cell_type, bed_paths in bed_files_by_celltype.items():
        if cell_type not in celltype_to_task_idx:
            raise KeyError(f"cell_type '{cell_type}' missing from celltype_to_task_idx")
        task_idx = int(celltype_to_task_idx[cell_type])

        for bed_path in bed_paths:
            label_df = load_label_intervals(bed_path)
            if label_df.empty:
                continue

            # Only consider chromosomes present in this BED
            for chrom in sorted(label_df["Chromosome"].unique()):
                tdf = transcript_df[transcript_df["Chromosome"] == chrom]
                if tdf.empty:
                    continue

                ldf_chr = label_df[label_df["Chromosome"] == chrom]

                for _chrom, w_start, w_end in iter_windows(tdf, window_bp, stride_bp):
                    seq = fasta.fetch(_chrom, w_start, w_end)
                    if "N" * window_bp == seq:
                        continue

                    if len(seq) != window_bp:
                        continue

                    y = np.zeros((n_bins, n_tasks), dtype=np.float32)

                    targ, counts = rasterize_window(
                        w_start, w_end, bin_size, ldf_chr, n_bins
                    )
                    y[:, task_idx] = targ

                    # Correct mask: 1 where there was *any* label support (transcript space / overlaps)
                    mask = (counts > 0).astype(np.float32)

                    examples.append({
                        "sequence": one_hot_encode(seq),
                        "targets": y,
                        "mask": mask,
                        "cell_type": cell_type,
                        "bed_path": str(bed_path),
                        "chrom": _chrom,
                        "window_start": int(w_start),
                        "window_end": int(w_end),
                    })

                    emitted += 1
                    if (max_examples is not None) and (emitted >= max_examples):
                        if examples:
                            save_shard(examples, out_dir, split, shard_id)
                        return

                    if len(examples) >= shard_size:
                        save_shard(examples, out_dir, split, shard_id)
                        shard_id += 1
                        examples = []

    if examples:
        save_shard(examples, out_dir, split, shard_id)


# ----------------------------
# CLI
# ----------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description="Convert per-sample BED soft labels into windowed (sequence, targets) shards."
    )
    ap.add_argument("--fasta", required=True, help="Reference genome FASTA (indexed).")
    ap.add_argument("--transcripts", required=True, help="TSV with Chromosome/Start/End transcript spans (0-based, half-open).")
    ap.add_argument("--out-dir", required=True, help="Output dataset directory.")
    ap.add_argument("--split", default="train", choices=["train", "val", "test"], help="Dataset split name.")
    ap.add_argument("--window-bp", type=int, default=131072)
    ap.add_argument("--bin-size", type=int, default=32)
    ap.add_argument("--stride-bp", type=int, default=65536)
    ap.add_argument("--shard-size", type=int, default=128)
    ap.add_argument("--max-examples", type=int, default=None, help="Cap examples for smoke tests.")

    # Repeated args: --cell-type B_Cell --bed path1 --bed path2 ...
    ap.add_argument(
        "--cell-type", action="append", default=[],
        help="Cell type name (repeat). Must correspond to one output head."
    )
    ap.add_argument(
        "--bed", action="append", default=[],
        help="BED path for the *most recent* --cell-type (repeat)."
    )
    return ap.parse_args()


def build_bed_map(cell_types: List[str], beds: List[str]) -> Dict[str, List[str]]:
    """
    Simple parser for paired --cell-type / --bed usage.
    Usage pattern:
      --cell-type B_Cell --bed a.bed --bed b.bed --cell-type T_Cell --bed c.bed
    """
    if not cell_types:
        raise ValueError("Provide at least one --cell-type")

    bed_map: Dict[str, List[str]] = {ct: [] for ct in cell_types}
    current = cell_types[0]
    ct_idx = 0

    # We cannot infer grouping from argparse alone; enforce a simple convention:
    # whenever --cell-type repeats, it advances; all --bed following belong to current until next cell-type.
    # Implement by scanning argv would be cleaner; here we require users to pass bed lists per ct
    # via comma-separated values as a safer interface.
    #
    # Since you’re operating in controlled scripts, easiest is to pass a manifest instead.
    raise NotImplementedError(
        "For reliability, use a manifest TSV. See below in response for recommended usage."
    )


def main():
    args = parse_args()

    fasta = FastaReader(args.fasta)
    transcript_df = load_transcripts_tsv(args.transcripts)

    manifest_path = Path(args.out_dir) / "manifest.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Expected manifest TSV at {manifest_path}"
        )

    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)

    required = {"cell_type", "bed_path"}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f"Manifest missing columns: {sorted(missing)}")

    bed_files_by_celltype: Dict[str, List[str]] = {}
    for ct, sub in manifest.groupby("cell_type"):
        bed_files_by_celltype[ct] = sub["bed_path"].tolist()

    celltypes = sorted(bed_files_by_celltype.keys())
    celltype_to_task_idx = {ct: i for i, ct in enumerate(celltypes)}

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
    )

if __name__ == "__main__":
    main()