import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import pyranges as pr

from generate_splicing_labels import (
    load_rmats_se,
    load_rmats_ri,
    load_rmats_alt_site,
    load_stringtie_gtf,
    _disjoint_labels_with_priority,   # REQUIRED
)

DEFAULTS = {
    "rmats_se": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_only/SE.MATS.JC.txt",
    "rmats_ri": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_only/RI.MATS.JC.txt",
    "rmats_a3ss": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_only/A3SS.MATS.JC.txt",
    "rmats_a5ss": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_only/A5SS.MATS.JC.txt",
    "gtf": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/compute_psi/stringtie_out/dominant_isoforms/B_Cell_dominant.gtf",
    "min_coverage": 20,
    "out_prefix": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/compute_psi/soft_labels/B_cell",
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Emit per-sample soft-label BEDs from rMATS output."
    )
    parser.add_argument("--rmats-se", default=DEFAULTS["rmats_se"])
    parser.add_argument("--rmats-ri", default=DEFAULTS["rmats_ri"])
    parser.add_argument("--rmats-a3ss", default=DEFAULTS["rmats_a3ss"])
    parser.add_argument("--rmats-a5ss", default=DEFAULTS["rmats_a5ss"])
    parser.add_argument("--gtf", default=DEFAULTS["gtf"])
    parser.add_argument("--min-coverage", type=int, default=DEFAULTS["min_coverage"])
    parser.add_argument("--out-prefix", default=DEFAULTS["out_prefix"])
    return parser.parse_args()

def get_psi_sample_columns(df):
    return sorted(c for c in df.columns if c.startswith("psi_"))

def build_shared_segments_df(stringtie_gr, rmats_all_bounds_df):
    """
    Build disjoint segments once, including exon/intron base label and Strand.
    Returns: Chromosome, Start, End, Strand, base_exon
    """
    gdf = stringtie_gr.df

    exons_df = gdf[gdf.Feature == "exon"][["Chromosome","Start","End","Strand"]].copy()
    transcripts_df = gdf[gdf.Feature == "transcript"][["Chromosome","Start","End","Strand"]].copy()

    if transcripts_df.empty and ("transcript_id" in gdf.columns):
        spans = (
            gdf[gdf.Feature == "exon"]
            .groupby(["Chromosome","Strand","transcript_id"], as_index=False)
            .agg(Start=("Start","min"), End=("End","max"))
        )
        transcripts_df = spans[["Chromosome","Start","End","Strand"]].copy()

    bounds = rmats_all_bounds_df[["Chromosome","Start","End","Strand"]].copy()
    bounds["psi"] = np.nan

    seg01 = _disjoint_labels_with_priority(
        transcripts_df=transcripts_df,
        exons_df=exons_df,
        psi_df=bounds
    )

    seg01 = seg01.rename(columns={"p_exon": "base_exon"})
    return seg01[["Chromosome","Start","End","Strand","base_exon"]]

def emit_beds_per_sample(
    stringtie_gr,
    rmats_se,
    rmats_ri,
    rmats_a3ss,
    rmats_a5ss,
    out_prefix
):
    rmats_all = pr.concat([rmats_se, rmats_ri, rmats_a3ss, rmats_a5ss])
    segments = build_shared_segments_df(stringtie_gr, rmats_all.df)

    assert segments[["Chromosome","Start","End","Strand"]].duplicated().sum() == 0

    psi_cols = get_psi_sample_columns(rmats_all.df)
    if not psi_cols:
        raise ValueError("No per-sample PSI columns found")

    print(f"[INFO] Found {len(psi_cols)} samples")

    for psi_col in psi_cols:
        print(f"[INFO] Emitting BED for {psi_col}")

        psi_events = (
            rmats_all.df[["Chromosome","Start","End","Strand", psi_col]]
            .rename(columns={psi_col: "psi"})
            .dropna(subset=["psi"])
        )

        bed = _disjoint_labels_with_priority(
            psi_df=psi_events,
            segments_df=segments
        )

        out_path = Path(f"{out_prefix}_{psi_col}.bed")
        out_path.parent.mkdir(parents=True, exist_ok=True)

        bed[["Chromosome","Start","End","p_exon"]].to_csv(
            out_path, sep="\t", header=False, index=False
        )

def main():
    args = parse_args()

    rmats_se = load_rmats_se(args.rmats_se, args.min_coverage)
    rmats_ri = load_rmats_ri(args.rmats_ri, args.min_coverage)
    rmats_a3ss = load_rmats_alt_site(args.rmats_a3ss, args.min_coverage)
    rmats_a5ss = load_rmats_alt_site(args.rmats_a5ss, args.min_coverage)

    stringtie_gr = load_stringtie_gtf(args.gtf)

    emit_beds_per_sample(
        stringtie_gr,
        rmats_se,
        rmats_ri,
        rmats_a3ss,
        rmats_a5ss,
        out_prefix=args.out_prefix
    )

if __name__ == "__main__":
    main()
