import pandas as pd
import numpy as np
import pyranges as pr

# --- Robust Helpers (Based on your fix) ---

def _normalize_strand(df_like):
    """
    Standardizes the Strand column to be a strict Category (+, -, .).
    This prevents PyRanges merge/concat failures.
    """
    if df_like is None:
        return None
    
    # Extract DataFrame
    if isinstance(df_like, pr.PyRanges):
        df = df_like.df.copy()
    else:
        df = df_like.copy()
        
    # Standardize Strand
    if 'Strand' not in df.columns:
        df['Strand'] = '.'
    
    # Force Categorical with explicit categories to align with PyRanges expectations
    df['Strand'] = pd.Categorical(df['Strand'].fillna('.'), categories=['+', '-', '.'])
    
    return pr.PyRanges(df)

def _ensure_pr(obj):
    """
    Guarantees the input is a valid, non-empty, normalized PyRanges object.
    Handles None, empty lists/DFs, and raw DataFrames.
    """
    # Handle None or Empty
    if obj is None:
        return pr.PyRanges()
    if hasattr(obj, '__len__') and len(obj) == 0:
        return pr.PyRanges()
        
    # Handle Existing PyRanges
    if isinstance(obj, pr.PyRanges):
        return _normalize_strand(obj)
        
    # Handle DataFrame
    if isinstance(obj, pd.DataFrame):
        required = {'Chromosome', 'Start', 'End'}
        missing = required - set(obj.columns)
        if missing:
            raise ValueError(f"Missing columns for PyRanges conversion: {missing}")
            
        obj = obj.copy()
        # Ensure integer coordinates
        obj['Start'] = obj['Start'].astype(int, copy=False)
        obj['End'] = obj['End'].astype(int, copy=False)
        return _normalize_strand(obj)
        
    raise TypeError(f"Unsupported type {type(obj)} for PyRanges conversion")

def _aggregate_exact_duplicates(df):
    """
    Resolves EXACT coordinate duplicates by averaging PSI and summing coverage.
    Includes explicit index reset to fix the 'Length of values' pandas crash.
    """
    if df.empty: return df
    
    # CRITICAL FIX: Reset index to ensure 0..N index.
    # PyRanges operations often leave fragmented/duplicate indices that crash groupby.
    df = df.reset_index(drop=True)
    
    agg_rules = {'psi': 'mean', 'coverage': 'sum'}
    
    # Use observed=True for categorical Strand efficiency
    # as_index=False returns a clean DataFrame directly
    return df.groupby(['Chromosome', 'Start', 'End', 'Strand'], observed=True, as_index=False).agg(agg_rules)

def _disjoint_labels_with_priority(transcripts_df=None, exons_df=None, psi_df=None, segments_df=None):
    """
    Two modes:

    1) Build + label (backward compatible):
        _disjoint_labels_with_priority(transcripts_df, exons_df, psi_df)
       Returns disjoint segments with p_exon assigned using priority PSI > exon > intron.

    2) Apply PSI onto PRECOMPUTED disjoint segments (new, for per-sample BEDs):
        _disjoint_labels_with_priority(segments_df=segments_df, psi_df=psi_df)
       segments_df must be disjoint and contain:
         Chromosome, Start, End, Strand, base_exon  (base_exon in {0.0, 1.0})
       Returns disjoint segments with p_exon assigned using priority PSI > base_exon.
    """
    def _prep(df, need_psi=False, need_base=False):
        if df is None or (hasattr(df, "empty") and df.empty):
            cols = ["Chromosome", "Start", "End", "Strand"]
            if need_psi:
                cols += ["psi"]
            if need_base:
                cols += ["base_exon"]
            return pd.DataFrame(columns=cols)

        df = df.copy()
        if "Strand" not in df.columns:
            df["Strand"] = "."
        df["Start"] = df["Start"].astype(int, copy=False)
        df["End"] = df["End"].astype(int, copy=False)

        keep = ["Chromosome", "Start", "End", "Strand"]
        if need_psi:
            keep += ["psi"]
        if need_base:
            keep += ["base_exon"]
        missing = [c for c in keep if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        return df[keep]

    # ---------------------------
    # MODE 2: Apply PSI to segments
    # ---------------------------
    if segments_df is not None:
        S = _prep(segments_df, need_base=True)
        P = _prep(psi_df, need_psi=True)

        if "Strand" not in S.columns:
            raise ValueError("segments_df must include Strand column")

        if S.empty:
            return pd.DataFrame(columns=["Chromosome","Start","End","p_exon","Strand"])


        # Ensure disjoint segments are sorted
        S = S.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
        P = P.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True) if not P.empty else P

        out_rows = []

        # Two-pointer overlap join per (chrom,strand)
        keys = S[["Chromosome", "Strand"]].drop_duplicates().itertuples(index=False, name=None)
        for chrom, strand in keys:
            seg = S[(S.Chromosome == chrom) & (S.Strand == strand)]
            psi = P[(P.Chromosome == chrom) & (P.Strand == strand)] if not P.empty else P


            seg_starts = seg["Start"].to_numpy()
            seg_ends = seg["End"].to_numpy()
            seg_base = seg["base_exon"].to_numpy()

            psi_starts = psi["Start"].to_numpy() if psi is not None and not psi.empty else np.array([], dtype=int)
            psi_ends = psi["End"].to_numpy() if psi is not None and not psi.empty else np.array([], dtype=int)
            psi_vals = psi["psi"].to_numpy() if psi is not None and not psi.empty else np.array([], dtype=float)

            j = 0  # pointer over psi events
            for i in range(len(seg)):
                s0 = int(seg_starts[i])
                s1 = int(seg_ends[i])

                # advance psi pointer while psi ends before segment starts
                while j < len(psi_ends) and int(psi_ends[j]) <= s0:
                    j += 1

                # collect all overlapping psi values
                k = j
                vals = []
                while k < len(psi_starts) and int(psi_starts[k]) < s1:
                    if int(psi_ends[k]) > s0:
                        vals.append(float(psi_vals[k]))
                    k += 1

                if vals:
                    p_exon = float(np.mean(vals))
                else:
                    p_exon = float(seg_base[i])

                out_rows.append((chrom, s0, s1, p_exon, strand))

        out = pd.DataFrame(out_rows, columns=["Chromosome", "Start", "End", "p_exon", "Strand"])
        return out

    # ---------------------------
    # MODE 1: Build disjoint segments + label (your existing behavior)
    # ---------------------------
    T = _prep(transcripts_df)
    E = _prep(exons_df)
    P = _prep(psi_df, need_psi=True)

    if T.empty:
        return pd.DataFrame(columns=["Chromosome", "Start", "End", "p_exon", "Strand"])

    out_rows = []

    keys = T[["Chromosome", "Strand"]].drop_duplicates().itertuples(index=False, name=None)
    for chrom, strand in keys:
        t = T[(T.Chromosome == chrom) & (T.Strand == strand)]
        e = E[(E.Chromosome == chrom) & (E.Strand == strand)]
        p = P[(P.Chromosome == chrom) & (P.Strand == strand)]

        boundaries = np.concatenate([
            t["Start"].to_numpy(), t["End"].to_numpy(),
            e["Start"].to_numpy() if not e.empty else np.array([], dtype=int),
            e["End"].to_numpy() if not e.empty else np.array([], dtype=int),
            p["Start"].to_numpy() if not p.empty else np.array([], dtype=int),
            p["End"].to_numpy() if not p.empty else np.array([], dtype=int),
        ])
        if boundaries.size == 0:
            continue
        boundaries = np.unique(boundaries)
        boundaries.sort()

        t_add, t_rem = {}, {}
        for s, en in zip(t["Start"].to_numpy(), t["End"].to_numpy()):
            t_add.setdefault(int(s), 0); t_add[int(s)] += 1
            t_rem.setdefault(int(en), 0); t_rem[int(en)] += 1

        e_add, e_rem = {}, {}
        if not e.empty:
            for s, en in zip(e["Start"].to_numpy(), e["End"].to_numpy()):
                e_add.setdefault(int(s), 0); e_add[int(s)] += 1
                e_rem.setdefault(int(en), 0); e_rem[int(en)] += 1

        p_add, p_rem = {}, {}
        if not p.empty:
            p = p.reset_index(drop=True)
            for idx, (s, en, psi) in enumerate(zip(p["Start"], p["End"], p["psi"])):
                s = int(s); en = int(en)
                p_add.setdefault(s, []).append((idx, float(psi)))
                p_rem.setdefault(en, []).append((idx, float(psi)))

        active_t = 0
        active_e = 0
        active_p_ids = {}
        active_p_sum = 0.0

        def _p_mean():
            n = len(active_p_ids)
            return (active_p_sum / n) if n else None

        prev_label = None
        prev_start = None

        for i in range(len(boundaries) - 1):
            x = int(boundaries[i])
            y = int(boundaries[i + 1])
            if x == y:
                continue

            if x in t_add: active_t += t_add[x]
            if x in t_rem: active_t -= t_rem[x]
            if x in e_add: active_e += e_add[x]
            if x in e_rem: active_e -= e_rem[x]
            if x in p_add:
                for pid, psi in p_add[x]:
                    if pid not in active_p_ids:
                        active_p_ids[pid] = psi
                        active_p_sum += psi
            if x in p_rem:
                for pid, _psi in p_rem[x]:
                    if pid in active_p_ids:
                        active_p_sum -= active_p_ids.pop(pid)

            if active_t <= 0:
                if prev_start is not None:
                    out_rows.append((chrom, prev_start, x, prev_label, strand))
                    prev_start, prev_label = None, None
                continue

            pm = _p_mean()
            label = pm if pm is not None else (1.0 if active_e > 0 else 0.0)

            if prev_start is None:
                prev_start = x
                prev_label = label
            elif label != prev_label:
                out_rows.append((chrom, prev_start, x, prev_label, strand))
                prev_start = x
                prev_label = label

        last_boundary = int(boundaries[-1])
        if prev_start is not None and prev_start < last_boundary:
            out_rows.append((chrom, prev_start, last_boundary, prev_label, strand))

    out = pd.DataFrame(out_rows, columns=["Chromosome", "Start", "End", "p_exon", "Strand"])
    return out[["Chromosome", "Start", "End", "p_exon", "Strand"]]


def _load_rmats_common(
    rmats_file,
    min_coverage,
    mode="per-sample",      # "per-sample" or "mean"
):
    """
    Load rMATS output while preserving per-sample PSI.

    mode:
      - "per-sample": keep psi_0, psi_1, ... columns
      - "mean": additionally compute psi_mean (averaged across non-NA samples)

    Returns a pandas DataFrame (NOT PyRanges yet).
    """
    try:
        df = pd.read_csv(rmats_file, sep="\t")
    except FileNotFoundError:
        print(f"[WARN] File not found: {rmats_file}")
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        print(f"[WARN] File empty: {rmats_file}")
        return pd.DataFrame()
    except pd.errors.ParserError as exc:
        print(f"[WARN] Failed to parse {rmats_file}: {exc}")
        return pd.DataFrame()

    # -----------------------------
    # Parse PSI per sample
    # -----------------------------
    def parse_float_list(values):
        if pd.isna(values):
            return []
        tokens = [t for t in str(values).split(",") if t not in {"", "NA"}]
        out = []
        for t in tokens:
            try:
                out.append(float(t))
            except ValueError:
                continue
        return out

    psi_lists = df["IncLevel1"].apply(parse_float_list)

    max_n = psi_lists.map(len).max()
    if pd.isna(max_n) or max_n == 0:
        return pd.DataFrame()

    psi_cols = [f"psi_{i}" for i in range(max_n)]

    psi_df = pd.DataFrame(
        psi_lists.tolist(),
        columns=psi_cols,
        index=df.index
    )

    df = pd.concat([df, psi_df], axis=1)

    if mode == "mean":
        df["psi_mean"] = psi_df.mean(axis=1, skipna=True)

    # -----------------------------
    # Parse coverage (unchanged)
    # -----------------------------
    count_cols = [c for c in df.columns if c.startswith(("IJC_SAMPLE", "SJC_SAMPLE"))]

    def parse_int_sum(values):
        if pd.isna(values):
            return 0
        total = 0
        for token in str(values).split(","):
            try:
                total += int(float(token))
            except ValueError:
                continue
        return total

    for col in count_cols:
        df[col] = df[col].apply(parse_int_sum)

    df["coverage"] = df[count_cols].sum(axis=1)

    # -----------------------------
    # Final filtering
    # -----------------------------
    df = df[df["coverage"] >= min_coverage]

    # Drop rows where *all* psi_i are NA
    psi_non_na = psi_df.notna().any(axis=1)
    df = df.loc[psi_non_na]
    df = df.reset_index(drop=True)

    return df

# --- Loaders ---

def load_rmats_se(rmats_file, min_coverage=20, mode="per-sample"):
    df = _load_rmats_common(rmats_file, min_coverage, mode=mode)
    if df.empty:
        return pr.PyRanges()

    df = df.rename(columns={
        'chr': 'Chromosome',
        'exonStart_0base': 'Start',
        'exonEnd': 'End',
        'strand': 'Strand'
    })
    df[['Start', 'End']] = df[['Start', 'End']].astype(int)
    df = df[df['Start'] < df['End']]

    return _normalize_strand(df)

def load_rmats_ri(rmats_file, min_coverage=20, mode="per-sample"):
    df = _load_rmats_common(rmats_file, min_coverage, mode=mode)
    if df.empty:
        return pr.PyRanges()

    df["Start"] = df["upstreamEE"]
    df["End"] = df["downstreamES"]
    df = df.rename(columns={
        "chr": "Chromosome",
        "strand": "Strand"
    })
    df[["Start", "End"]] = df[["Start", "End"]].astype(int)
    df = df[df["Start"] < df["End"]]

    return _normalize_strand(df)

def load_rmats_alt_site(rmats_file, min_coverage=20, mode="per-sample"):
    df = _load_rmats_common(rmats_file, min_coverage, mode=mode)
    if df.empty:
        return pr.PyRanges()

    # Long exon
    long_df = df.rename(columns={
        'chr': 'Chromosome',
        'longExonStart_0base': 'Start',
        'longExonEnd': 'End',
        'strand': 'Strand'
    })[['Chromosome','Start','End','Strand']]

    # Short exon
    short_df = df.rename(columns={
        'chr': 'Chromosome',
        'shortES': 'Start',
        'shortEE': 'End',
        'strand': 'Strand'
    })[['Chromosome','Start','End','Strand']]

    pr_long = _normalize_strand(long_df)
    pr_short = _normalize_strand(short_df)

    extension = pr_long.subtract(pr_short)

    if len(extension) == 0:
        return pr.PyRanges()

    return _normalize_strand(extension.df)

def load_stringtie_gtf(gtf_file):
    return pr.read_gtf(gtf_file)

# --- Generators ---

def generate_hard_labels(stringtie_gr, rmats_se, rmats_ri=None, rmats_a3ss=None, rmats_a5ss=None, psi_low=0.2, psi_high=0.8):
    # Ensure all inputs are robust PyRanges
    stringtie_gr = _ensure_pr(stringtie_gr)
    rmats_se = _ensure_pr(rmats_se)
    rmats_ri = _ensure_pr(rmats_ri)
    rmats_a3ss = _ensure_pr(rmats_a3ss)
    rmats_a5ss = _ensure_pr(rmats_a5ss)

    st_exons = stringtie_gr[stringtie_gr.Feature == 'exon']
    st_transcripts = stringtie_gr[stringtie_gr.Feature == 'transcript']
    
    if len(st_transcripts) == 0:
        st_transcripts = st_exons.merge(by="transcript_id")

    # 1. Correct Exons (SE)
    false_exons = rmats_se[rmats_se.psi < psi_low]
    corrected_exons = st_exons.subtract(false_exons) if len(false_exons) else st_exons

    missed_exons = rmats_se[rmats_se.psi > psi_high]
    if len(missed_exons):
        corrected_exons = pr.concat([corrected_exons, missed_exons]).merge()

    # 2. Retained Introns (RI)
    if len(rmats_ri):
        ri_high = rmats_ri[rmats_ri.psi > psi_high]
        if len(ri_high):
            ri_high = _normalize_strand(ri_high)
            corrected_exons = pr.concat([corrected_exons, ri_high]).merge()

    # 3. Extensions (A3SS/A5SS)
    extensions_list = []
    if len(rmats_a3ss): extensions_list.append(rmats_a3ss)
    if len(rmats_a5ss): extensions_list.append(rmats_a5ss)
    
    extensions = pr.concat(extensions_list) if extensions_list else pr.PyRanges()

    if len(extensions):
        ext_high = extensions[extensions.psi > psi_high]
        if len(ext_high):
            ext_high = _normalize_strand(ext_high)
            corrected_exons = pr.concat([corrected_exons, ext_high]).merge()
            
        ext_low = extensions[extensions.psi < psi_low]
        if len(ext_low):
            ext_low = _normalize_strand(ext_low)
            corrected_exons = corrected_exons.subtract(ext_low)

    # 4. Define Introns
    introns = st_transcripts.subtract(corrected_exons)

    corrected_exons = _normalize_strand(corrected_exons)
    introns = _normalize_strand(introns)

    corrected_exons.label = 2
    introns.label = 1
    
    return pr.concat([corrected_exons, introns]).sort().df

def generate_soft_labels(stringtie_gr, rmats_se, rmats_ri=None, rmats_a3ss=None, rmats_a5ss=None):
    # Ensure all inputs are robust PyRanges
    stringtie_gr = _ensure_pr(stringtie_gr)
    rmats_se = _ensure_pr(rmats_se)
    rmats_ri = _ensure_pr(rmats_ri)
    rmats_a3ss = _ensure_pr(rmats_a3ss)
    rmats_a5ss = _ensure_pr(rmats_a5ss)

    # PSI regions: concat + exact-duplicate aggregation (your current behavior)
    vars_list = [gr for gr in (rmats_se, rmats_ri, rmats_a3ss, rmats_a5ss) if len(gr) > 0]
    vars_pr = pr.concat(vars_list) if vars_list else pr.PyRanges()
    psi_df = pd.DataFrame(columns=['Chromosome','Start','End','Strand','psi'])
    if len(vars_pr):
        vars_df = _aggregate_exact_duplicates(vars_pr.df)  # averages PSI for exact duplicates
        psi_df = vars_df[['Chromosome','Start','End','Strand','psi']].copy()

    # Transcript and exon space from the dominant-isoform GTF
    gdf = stringtie_gr.df
    exons_df = gdf[gdf.Feature == 'exon'][['Chromosome','Start','End','Strand']].copy()
    transcripts_df = gdf[gdf.Feature == 'transcript'][['Chromosome','Start','End','Strand']].copy()
    if transcripts_df.empty:
        # fallback: build transcript spans from exons
        # (matches your prior behavior of merging by transcript_id, but avoids PyRanges ops here)
        if 'transcript_id' in gdf.columns and not exons_df.empty:
            spans = gdf[gdf.Feature == 'exon'].groupby(
                ['Chromosome','Strand','transcript_id'], as_index=False
            ).agg(Start=('Start','min'), End=('End','max'))
            transcripts_df = spans[['Chromosome','Start','End','Strand']].copy()
        else:
            transcripts_df = exons_df.copy()

    # Build disjoint labeled segments with PSI priority
    out = _disjoint_labels_with_priority(transcripts_df, exons_df, psi_df)
    return out
