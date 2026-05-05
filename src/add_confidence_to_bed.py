import argparse
from pathlib import Path
import logging

import pandas as pd
import numpy as np
import pyranges as pr

# We still import these for reference or if we fall back, but we will override the loaders
# to include confidence parsing.
from generate_splicing_labels import (
    load_stringtie_gtf,
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
        description="Emit per-sample soft-label BEDs from rMATS output with Confidence."
    )
    parser.add_argument("--rmats-se", default=DEFAULTS["rmats_se"])
    parser.add_argument("--rmats-ri", default=DEFAULTS["rmats_ri"])
    parser.add_argument("--rmats-a3ss", default=DEFAULTS["rmats_a3ss"])
    parser.add_argument("--rmats-a5ss", default=DEFAULTS["rmats_a5ss"])
    parser.add_argument("--gtf", default=DEFAULTS["gtf"])
    parser.add_argument("--min-coverage", type=int, default=DEFAULTS["min_coverage"])
    parser.add_argument("--out-prefix", default=DEFAULTS["out_prefix"])
    return parser.parse_args()

# =============================================================================
# Custom Loaders to Parse IJC/SJC per replicate for Confidence Calculation
# =============================================================================

def _normalize_strand(df_like):
    if df_like is None:
        return None
    if isinstance(df_like, pr.PyRanges):
        df = df_like.df.copy()
    else:
        df = df_like.copy()
    if 'Strand' not in df.columns:
        df['Strand'] = '.'
    df['Strand'] = pd.Categorical(df['Strand'].fillna('.'), categories=['+', '-', '.'])
    return pr.PyRanges(df)

def _load_rmats_common_custom(rmats_file, min_coverage=20, mode="per-sample"):
    try:
        df = pd.read_csv(rmats_file, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        print(f"[WARN] Failed to load {rmats_file}: {e}")
        return pd.DataFrame()

    def parse_float_list(values):
        if pd.isna(values): return []
        tokens = str(values).split(",")
        return [float(t) if t not in {"", "NA"} else np.nan for t in tokens]

    def parse_int_list(values):
        if pd.isna(values): return []
        tokens = str(values).split(",")
        return [int(float(t)) if t not in {"", "NA"} else 0 for t in tokens]

    # 1. Parse PSI (IncLevel1) -> psi_0, psi_1...
    psi_lists = df["IncLevel1"].apply(parse_float_list)
    max_n = psi_lists.map(len).max()
    if pd.isna(max_n) or max_n == 0:
        return pd.DataFrame()
    
    psi_cols = [f"psi_{i}" for i in range(max_n)]
    psi_df = pd.DataFrame(psi_lists.tolist(), columns=psi_cols, index=df.index)
    
    # 2. Parse IJC (IJC_SAMPLE_1) -> ijc_0, ijc_1...
    ijc_lists = df["IJC_SAMPLE_1"].apply(parse_int_list)
    ijc_cols = [f"ijc_{i}" for i in range(max_n)]
    ijc_df = pd.DataFrame(ijc_lists.tolist(), columns=ijc_cols, index=df.index)

    # 3. Parse SJC (SJC_SAMPLE_1) -> sjc_0, sjc_1...
    sjc_lists = df["SJC_SAMPLE_1"].apply(parse_int_list)
    sjc_cols = [f"sjc_{i}" for i in range(max_n)]
    sjc_df = pd.DataFrame(sjc_lists.tolist(), columns=sjc_cols, index=df.index)

    # Combine
    df = pd.concat([df, psi_df, ijc_df, sjc_df], axis=1)

    # 4. Coverage Filtering
    # To match original logic, we filter based on TOTAL coverage in SAMPLE 1
    # Original used: df[count_cols].sum(axis=1) where count_cols included SAMPLE_1 and SAMPLE_2(if exists)
    # We will compute coverage for filtering based on IJC_1 + SJC_1 (Total reads for this group)
    # Note: If min_coverage logic previously included Sample 2, we might be stricter here if we only look at Sample 1.
    # But usually we care if *this* sample has coverage.
    
    # Correctly sum up all IJC/SJC columns
    total_cov = ijc_df.fillna(0).sum(axis=1) + sjc_df.fillna(0).sum(axis=1)
    
    df = df[total_cov >= min_coverage]

    # Drop rows where all PSIs are NA
    psi_non_na = psi_df.notna().any(axis=1)
    df = df.loc[psi_non_na].reset_index(drop=True)
    
    return df

def load_rmats_se_custom(rmats_file, min_coverage=20):
    df = _load_rmats_common_custom(rmats_file, min_coverage)
    if df.empty: return pr.PyRanges()
    df = df.rename(columns={'chr': 'Chromosome', 'exonStart_0base': 'Start', 'exonEnd': 'End', 'strand': 'Strand'})
    df[['Start', 'End']] = df[['Start', 'End']].astype(int)
    df = df[df['Start'] < df['End']]
    return _normalize_strand(df)

def load_rmats_ri_custom(rmats_file, min_coverage=20):
    df = _load_rmats_common_custom(rmats_file, min_coverage)
    if df.empty: return pr.PyRanges()
    df["Start"] = df["upstreamEE"]
    df["End"] = df["downstreamES"]
    df = df.rename(columns={"chr": "Chromosome", "strand": "Strand"})
    df[["Start", "End"]] = df[["Start", "End"]].astype(int)
    df = df[df["Start"] < df["End"]]
    return _normalize_strand(df)

def load_rmats_alt_site_custom(rmats_file, min_coverage=20):
    df = _load_rmats_common_custom(rmats_file, min_coverage)
    if df.empty: return pr.PyRanges()
    
    long_df = df.rename(columns={'chr': 'Chromosome', 'longExonStart_0base': 'Start', 'longExonEnd': 'End', 'strand': 'Strand'})
    short_df = df.rename(columns={'chr': 'Chromosome', 'shortES': 'Start', 'shortEE': 'End', 'strand': 'Strand'})
    
    # Select cols to keep (Coordinates + Payload)
    # Payload = psi_*, ijc_*, sjc_*
    payload_cols = [c for c in df.columns if c.startswith(("psi_", "ijc_", "sjc_"))]
    cols = ['Chromosome', 'Start', 'End', 'Strand'] + payload_cols
    
    pr_long = _normalize_strand(long_df[cols])
    pr_short = _normalize_strand(short_df[cols])
    
    # Subtracting short from long to get the 'extension'
    # Default PyRanges subtraction keeps metadata of the "self" (long)
    extension = pr_long.subtract(pr_short)
    
    if len(extension) == 0:
        return pr.PyRanges()
        
    return _normalize_strand(extension.df)


# =============================================================================
# Custom Disjoint Labeling with Confidence Aggregation
# =============================================================================

def _disjoint_labels_with_priority_custom(transcripts_df=None, exons_df=None, psi_df=None, segments_df=None):
    """
    Custom version that propagates 'confidence_value' alongside 'psi'.
    Expects psi_df to have ['psi', 'confidence_value'].
    """
    # ... Helper _prep ...
    def _prep(df, need_psi=False, need_base=False):
        if df is None or (hasattr(df, "empty") and df.empty):
            cols = ["Chromosome", "Start", "End", "Strand"]
            if need_psi: cols += ["psi", "confidence_value"]
            if need_base: cols += ["base_exon"]
            return pd.DataFrame(columns=cols)
        return df.copy()

    # Mode 2: Apply PSI + Confidence to Segments
    if segments_df is not None:
        T = pd.DataFrame(columns=["Chromosome", "Start", "End", "Strand"])
        E = pd.DataFrame(columns=["Chromosome", "Start", "End", "Strand"])
        P = _prep(psi_df, need_psi=True)
        
        # We treat segments as the "Transcripts" backbone for iterating
        # Actually, we want to intersect P onto Segments.
        # But to reuse the sweep-line logic, we can pass segments as T/E?
        # No, the original logic for Mode 2 uses segments as T?
        
        # Actually, checking original code (not visible, but inferred):
        # The logic below performs a union of boundaries.
        # If we pass segments_df as 'exons_df' effectively (or just a backbone),
        # we can attribute P over S.
        
        # Let's adapt the sweep line exactly for this purpose:
        # We iterate over Chrom/Strand.
        # Valid intervals are defined by segments_df.
        # PSI/Confidence comes from psi_df.
        
        # To strictly map PSI to Precomputed Segments without splitting them further:
        # We can use PyRanges join?
        # Or just use the sweep line which handles overlaps naturally.
        
        # Let's use the sweep line but ensure we output the segments 
        # enriched with the mean PSI/Confidence.
        
        S = segments_df.copy()
        P = P.dropna(subset=["psi"])
        
        if P.empty:
            S["p_exon"] = S["base_exon"]
            S["confidence_value"] = 0.0 # Default confidence for background
            return S
        
        # Use PyRanges to map efficiently? 
        # Sweep line is good for averaging overlaps.
        
        out_rows = []
        keys = S[["Chromosome", "Strand"]].drop_duplicates().itertuples(index=False, name=None)
        
        for chrom, strand in keys:
            # S_sub: The segments (targets)
            s_sub = S[(S.Chromosome == chrom) & (S.Strand == strand)]
            # P_sub: The data
            p_sub = P[(P.Chromosome == chrom) & (P.Strand == strand)]

            # We process this by just using PyRanges join/overlap if we trust segments are disjoint
            # However, multiple PSIs can overlap one segment. We want weighted average?
            # Or just average.
            
            # Using a simplified sweep or interval tree for P over S.
            # Since S is disjoint, for each segment in S, find overlapping Ps.
            
            # Using PyRanges for speed
            pr_s = pr.PyRanges(s_sub)
            pr_p = pr.PyRanges(p_sub)
            
            # Join P to S
            # S is "self", P is "other". 
            # We want to annotate S with P data.
            j = pr_s.join(pr_p, report_overlap=True)
            
            if j.empty:
                # No overlaps in this chrom/strand
                # Keep base exon, conf=0
                sub_res = s_sub.copy()
                sub_res["p_exon"] = sub_res["base_exon"]
                sub_res["confidence_value"] = 0.0
                out_rows.append(sub_res)
                continue
                
            dfj = j.df
            
            # Weighted average of PSI and Confidence by Overlap length
            dfj["w_psi"] = dfj["psi"] * dfj["Overlap"]
            dfj["w_conf"] = dfj["confidence_value"] * dfj["Overlap"]
            
            agg = dfj.groupby(["Chromosome", "Start", "End", "Strand"], observed=True).agg(
                sum_w_psi=("w_psi", "sum"),
                sum_w_conf=("w_conf", "sum"),
                total_overlap=("Overlap", "sum")
            ).reset_index()
            
            agg["p_exon"] = agg["sum_w_psi"] / agg["total_overlap"]
            agg["confidence_value"] = agg["sum_w_conf"] / agg["total_overlap"]
            
            # Merge back to original S to catch those with NO overlap
            res = pd.merge(
                s_sub, 
                agg[["Chromosome", "Start", "End", "Strand", "p_exon", "confidence_value"]],
                on=["Chromosome", "Start", "End", "Strand"],
                how="left"
            )
            
            # Fill NA
            res["p_exon"] = res["p_exon"].fillna(res["base_exon"])
            res["confidence_value"] = res["confidence_value"].fillna(0.0)
            
            out_rows.append(res)
            
        if not out_rows:
            return pd.DataFrame()
            
        return pd.concat(out_rows).reset_index(drop=True)

    # Mode 1 (Original Build) is not needed for this script if we only call with segments_df
    # But for build_shared_segments_df, we need Mode 1.
    # We can perform a dummy call to the original logical or implement a simplified logic 
    # since we only need Mode 1 for `build_shared_segments_df` which doesn't use confidence.
    pass      

# We need the original Mode 1 logic for 'build_shared_segments_df' to work
# Let's import the original for Mode 1 use, and use our custom one for Mode 2.
from generate_splicing_labels import _disjoint_labels_with_priority as _disjoint_labels_original


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

    # Use ORIGNAL logic for this structural step
    seg01 = _disjoint_labels_original(
        transcripts_df=transcripts_df,
        exons_df=exons_df,
        psi_df=bounds
    )

    seg01 = seg01.rename(columns={"p_exon": "base_exon"})
    return seg01[["Chromosome","Start","End","Strand","base_exon"]]


def get_psi_sample_columns(df):
    return sorted(c for c in df.columns if c.startswith("psi_"))

def emit_beds_per_sample(
    stringtie_gr,
    rmats_se,
    rmats_ri,
    rmats_a3ss,
    rmats_a5ss,
    out_prefix
):
    rmats_all = pr.concat([rmats_se, rmats_ri, rmats_a3ss, rmats_a5ss])

    if len(rmats_all) == 0:
        print("[WARN] No rMATS events found across all categories. Skipping.")
        return

    segments = build_shared_segments_df(stringtie_gr, rmats_all.df)

    assert segments[["Chromosome","Start","End","Strand"]].duplicated().sum() == 0

    psi_cols = get_psi_sample_columns(rmats_all.df)
    if not psi_cols:
        raise ValueError("No per-sample PSI columns found")

    print(f"[INFO] Found {len(psi_cols)} samples: {psi_cols}")

    for psi_col in psi_cols:
        # Extract the index ID (e.g. "psi_0" -> "0")
        idx = psi_col.split("_")[1]
        ijc_col = f"ijc_{idx}"
        sjc_col = f"sjc_{idx}"
        
        print(f"[INFO] Processing {psi_col} (Reads: {ijc_col}, {sjc_col})")

        # Select relevant columns
        cols = ["Chromosome","Start","End","Strand", psi_col]
        # Check if ijc/sjc cols exist
        if ijc_col in rmats_all.df.columns and sjc_col in rmats_all.df.columns:
            cols.extend([ijc_col, sjc_col])
        else:
            print(f"[WARN] Coverage columns {ijc_col}/{sjc_col} not found. Confidence set to 0.")
        
        psi_events = rmats_all.df[cols].copy()
        psi_events = psi_events.rename(columns={psi_col: "psi"})
        psi_events = psi_events.dropna(subset=["psi"])
        
        # Calculate Confidence
        # log(1 + IJC + SJC)
        if ijc_col in psi_events.columns:
            ijc_vals = psi_events[ijc_col].fillna(0)
            sjc_vals = psi_events[sjc_col].fillna(0)
            psi_events["confidence_value"] = np.log1p(ijc_vals + sjc_vals)
        else:
            psi_events["confidence_value"] = 0.0

        # Use CUSTOM logic for mapping to segments
        bed = _disjoint_labels_with_priority_custom(
            psi_df=psi_events,
            segments_df=segments
        )

        out_path = Path(f"{out_prefix}_{psi_col}.bed")
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Output additional confidence column
        bed[["Chromosome","Start","End","p_exon", "confidence_value"]].to_csv(
            out_path, sep="\t", header=False, index=False
        )
        print(f"       -> Wrote {out_path}")

def main():
    args = parse_args()

    # Use CUSTOM loaders
    print("[INFO] Loading rMATS files with custom confidence parser...")
    rmats_se = load_rmats_se_custom(args.rmats_se, args.min_coverage)
    rmats_ri = load_rmats_ri_custom(args.rmats_ri, args.min_coverage)
    rmats_a3ss = load_rmats_alt_site_custom(args.rmats_a3ss, args.min_coverage)
    rmats_a5ss = load_rmats_alt_site_custom(args.rmats_a5ss, args.min_coverage)

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
