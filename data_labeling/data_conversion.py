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
    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df = df.dropna(subset=["Start", "End"]).reset_index(drop=True)
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

    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df = df.dropna(subset=["Start", "End"]).reset_index(drop=True)

    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df = df[df["Start"] < df["End"]].reset_index(drop=True)
    return df


# ----------------------------
# Windowing + rasterization
# ----------------------------

# NOTE: Strand is intentionally ignored.
# Borzoi learns strand via RC augmentation.

def iter_windows_from_regions(
    regions_df: pd.DataFrame,
    window_bp: int,
    stride_bp: int,
):
    """
    Yield genomic windows from the UNION of transcript regions per chromosome.

    regions_df must contain Chromosome, Start, End (0-based, half-open).
    """
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
    Length-weighted mean p_exon per bin.

    Returns:
      targets[b] = mean p_exon
      counts[b]  = bp covered by *any* transcript segment
      events[b]  = bp covered by a PSI-defined event (non-NaN)
    """
    targets = np.zeros(n_bins, dtype=np.float32)
    counts = np.zeros(n_bins, dtype=np.float32)
    events = np.zeros(n_bins, dtype=np.float32)

    if label_df.empty:
        return targets, counts, events

    for _, row in label_df.iterrows():
        s = max(int(row.Start), window_start)
        e = min(int(row.End), window_end)
        if s >= e:
            continue

        b0 = (s - window_start) // bin_size
        b1 = (e - window_start + bin_size - 1) // bin_size

        p = row.p_exon
        has_event = not np.isnan(p)

        for b in range(int(b0), int(b1)):
            bin_s = window_start + b * bin_size
            bin_e = bin_s + bin_size
            ov = max(0, min(e, bin_e) - max(s, bin_s))
            if ov <= 0:
                continue

            counts[b] += ov

            if has_event:
                targets[b] += ov * float(p)
                events[b] += ov

    nz = events > 0
    targets[nz] /= events[nz]

    return targets, counts, events


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
):
    if window_bp % bin_size != 0:
        raise ValueError(f"window_bp ({window_bp}) must be divisible by bin_size ({bin_size})")

    n_bins = window_bp // bin_size
    n_tasks = len(celltype_to_task_idx)

    required = {"Chromosome", "Start", "End"}
    missing = required - set(transcript_df.columns)
    if missing:
        raise ValueError(f"Transcript DF missing columns: {sorted(missing)}")

    examples: List[dict] = []
    shard_id = 0
    emitted = 0

    # Use transcript UNION for windowing
    FLANK = 262_144  # half-window

    # regions_df = transcript_df[["Chromosome","Start","End"]].copy()

    regions_df = transcript_df.copy()
    regions_df["Start"] = (regions_df["Start"] - FLANK).clip(lower=0)
    regions_df["End"] = regions_df["End"] + FLANK
    # Expand transcript spans so PSI events can land in the central prediction region
    # after Borzoi cropping (±163,840 bp) within a 524,288 bp input window.
    regions_df = regions_df[["Chromosome","Start","End"]]

    windows = list(iter_windows_from_regions(regions_df, window_bp, stride_bp))
    print(f"Total windows to scan: {len(windows)}")

    for cell_type, bed_paths in bed_files_by_celltype.items():
        task_idx = int(celltype_to_task_idx[cell_type])

        for bed_path in bed_paths:
            label_df = load_label_intervals(bed_path)
            if label_df.empty:
                continue

            for chrom, w_start, w_end in windows:
                # Subset labels to current chrom
                ldf_chr = label_df[label_df["Chromosome"] == chrom]
                if ldf_chr.empty:
                    continue

                seq = fasta.fetch(chrom, w_start, w_end)
                if len(seq) != window_bp:
                    continue

                # Skip windows that are mostly Ns (>50%)
                if seq.upper().count("N") / window_bp > 0.5:
                    continue

                y = np.zeros((n_bins, n_tasks), dtype=np.float32)

                targ, counts, events = rasterize_window(
                    w_start, w_end, bin_size, ldf_chr, n_bins
                )
                y[:, task_idx] = targ

                if np.any(y[:, np.arange(n_tasks) != task_idx]):
                    raise RuntimeError("Nonzero values detected in non-target task columns")

                # Mask indicates bins with PSI-defined supervision.
                # Bins outside PSI events contribute zero loss.
                mask = (events > 0).astype(np.float32)

                has_signal = float(mask.sum() > 0)

                examples.append({
                    "sequence": one_hot_encode(seq),
                    "targets": y,
                    "mask": mask,
                    "cell_type": cell_type,
                    "bed_path": str(bed_path),
                    "chrom": chrom,
                    "window_start": int(w_start),
                    "window_end": int(w_end),
                    "has_signal": has_signal,
                })

                emitted += 1
                if max_examples is not None and emitted >= max_examples:
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
    ap.add_argument(
        "--manifest",
        default=None,
        help=(
            "Optional TSV manifest mapping cell_type to bed_path. "
            "If omitted, defaults to <out-dir>/manifest.tsv."
        ),
    )
    ap.add_argument(
        "--chroms",
        type=str,
        default=None,
        help="Comma-separated list of chromosomes to include (e.g. chr1,chr2)"
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

    if args.chroms is not None:
        keep = set(c.strip() for c in args.chroms.split(","))
        transcript_df = transcript_df[transcript_df["Chromosome"].isin(keep)]

    if args.manifest is not None:
        manifest_path = Path(args.manifest)
    else:
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

    # Manifest is expected to describe a single multitask dataset covering all cell types.
    bed_files_by_celltype: Dict[str, List[str]] = {}
    for ct, sub in manifest.groupby("cell_type", sort=False):
        if ct not in bed_files_by_celltype:
            bed_files_by_celltype[ct] = []
        bed_files_by_celltype[ct].extend(sub["bed_path"].tolist())

    celltypes = sorted(bed_files_by_celltype.keys())
    celltype_to_task_idx = {ct: i for i, ct in enumerate(celltypes)}

    print("Building multitask dataset with cell types in order:")
    for ct, idx in celltype_to_task_idx.items():
        print(f"  task[{idx}] = {ct}")

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