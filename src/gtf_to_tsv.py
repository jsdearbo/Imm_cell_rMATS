#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import pyranges as pr


def gtf_to_transcripts_tsv(
    gtf_path: str,
    out_tsv: str,
    require_transcript_feature: bool = False,
    union_out: str = None,
):
    """
    Convert a dominant-isoform GTF to a transcripts.tsv file.

    Output coordinates are 0-based, half-open.
    """

    gr = pr.read_gtf(gtf_path)
    df = gr.df.copy()

    # ----------------------------------------
    # Prefer transcript features if present
    # ----------------------------------------
    has_transcripts = (df.Feature == "transcript").any()

    if has_transcripts:
        tx = df[df.Feature == "transcript"].copy()
    else:
        if require_transcript_feature:
            raise RuntimeError(
                "GTF has no 'transcript' features and require_transcript_feature=True"
            )

        # ----------------------------------------
        # Fallback: build transcript spans from exons
        # ----------------------------------------
        if "transcript_id" not in df.columns:
            raise RuntimeError(
                "GTF has no transcript features and no transcript_id column"
            )

        exons = df[df.Feature == "exon"].copy()
        if exons.empty:
            raise RuntimeError("No exon features found in GTF")

        tx = (
            exons
            .groupby(
                ["Chromosome", "Strand", "transcript_id", "gene_id"],
                observed=True,
                as_index=False
            )
            .agg(
                Start=("Start", "min"),
                End=("End", "max")
            )
        )

    # ----------------------------------------
    # Ensure required columns
    # ----------------------------------------
    required = {"Chromosome", "Start", "End", "Strand"}
    missing = required - set(tx.columns)
    if missing:
        raise RuntimeError(f"Missing required columns: {sorted(missing)}")

    # ----------------------------------------
    # Coordinate sanity
    # StringTie GTF is 1-based inclusive
    # PyRanges converts to 0-based half-open already
    # ----------------------------------------
    tx["Start"] = tx["Start"].astype(int)
    tx["End"] = tx["End"].astype(int)

    # Drop zero-length or invalid transcripts
    tx = tx[tx["Start"] < tx["End"]].copy()

    # ----------------------------------------
    # Keep a clean, minimal schema
    # ----------------------------------------
    keep_cols = ["Chromosome", "Start", "End", "Strand"]
    for optional in ["transcript_id", "gene_id", "gene_name"]:
        if optional in tx.columns:
            keep_cols.append(optional)

    tx = tx[keep_cols].sort_values(
        ["Chromosome", "Start", "End", "Strand"]
    ).reset_index(drop=True)

    # ----------------------------------------
    # Write TSV
    # ----------------------------------------
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    tx.to_csv(out_tsv, sep="\t", index=False)

    print(f"[INFO] Wrote {len(tx)} transcripts to {out_tsv}")

    if union_out is not None:
        union_path = Path(union_out)
        union_path.parent.mkdir(parents=True, exist_ok=True)

        if union_path.exists():
            existing = pd.read_csv(union_path, sep="\t")
        else:
            existing = pd.DataFrame(columns=tx.columns)

        combined = (
            pd.concat([existing, tx], ignore_index=True)
            .drop_duplicates()
            .sort_values(tx.columns.tolist())
            .reset_index(drop=True)
        )

        combined.to_csv(union_path, sep="\t", index=False)
        print(f"[INFO] Updated union transcripts TSV at {union_path} ({len(combined)} unique rows)")


def parse_args():
    ap = argparse.ArgumentParser(
        description="Convert dominant-isoform GTF to transcripts.tsv (0-based, half-open)."
    )
    ap.add_argument("--gtf", required=True, help="Dominant isoform GTF")
    ap.add_argument("--out", required=True, help="Output transcripts.tsv")
    ap.add_argument(
        "--require-transcript-feature",
        action="store_true",
        help="Fail if GTF lacks transcript features (no exon fallback)",
    )
    ap.add_argument(
        "--union-out",
        default=None,
        help="Optional path to write/update the union transcripts TSV across invocations.",
    )
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    gtf_to_transcripts_tsv(
        gtf_path=args.gtf,
        out_tsv=args.out,
        require_transcript_feature=args.require_transcript_feature,
        union_out=args.union_out,
    )
