import pandas as pd
import numpy as np
import pyranges as pr

# --- Robust Helpers ---

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
    if obj is None:
        return pr.PyRanges()
    if hasattr(obj, '__len__') and len(obj) == 0:
        return pr.PyRanges()
        
    if isinstance(obj, pr.PyRanges):
        return _normalize_strand(obj)
        
    if isinstance(obj, pd.DataFrame):
        required = {'Chromosome', 'Start', 'End'}
        missing = required - set(obj.columns)
        if missing:
            raise ValueError(f"Missing columns for PyRanges conversion: {missing}")
            
        obj = obj.copy()
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
    df = df.reset_index(drop=True)
    
    agg_rules = {}
    if 'psi' in df.columns: agg_rules['psi'] = 'mean'
    if 'coverage' in df.columns: agg_rules['coverage'] = 'sum'
    
    # Use observed=True for categorical Strand efficiency
    if not agg_rules:
        return df.drop_duplicates(subset=['Chromosome', 'Start', 'End', 'Strand'])

    return df.groupby(['Chromosome', 'Start', 'End', 'Strand'], observed=True, as_index=False).agg(agg_rules)

def _resolve_overlaps_mean_psi(gr):
    """
    Resolves overlapping intervals by splitting them into disjoint 'atomic' segments
    and averaging the PSI values.
    """
    if len(gr) == 0: return gr
    
    df = gr.df.reset_index(drop=True)
    results = []
    
    # Process per Chromosome/Strand group
    for (chrom, strand), group in df.groupby(['Chromosome', 'Strand'], observed=True):
        if group.empty: continue
        
        # Identify Atomic Intervals
        boundaries = np.unique(np.concatenate([group.Start.values, group.End.values]))
        if len(boundaries) < 2: continue
        
        starts = boundaries[:-1]
        ends = boundaries[1:]
        valid_mask = starts < ends
        starts = starts[valid_mask]
        ends = ends[valid_mask]
        
        if len(starts) == 0: continue

        elem_df = pd.DataFrame({
            'Chromosome': chrom,
            'Start': starts,
            'End': ends,
            'Strand': strand
        })
        elem_pr = pr.PyRanges(elem_df)
        group_pr = pr.PyRanges(group)
        
        # Join to map values
        joined = elem_pr.join(group_pr)
        if joined is None or len(joined) == 0: continue
        
        # Aggregate PSI
        agg = joined.df.groupby(['Chromosome', 'Start', 'End', 'Strand'], observed=True)['psi'].mean().reset_index()
        results.append(agg)
        
    if not results:
        return pr.PyRanges()
        
    final_df = pd.concat(results)
    return _normalize_strand(pr.PyRanges(final_df))

def _load_rmats_common(rmats_file, min_coverage):
    try:
        df = pd.read_csv(rmats_file, sep='\t')
    except Exception as exc:
        print(f"[WARN] Failed to parse {rmats_file}: {exc}")
        return pd.DataFrame()

    def parse_float_list(values):
        tokens = [t for t in str(values).split(',') if t not in {'', 'NA'}]
        return [float(t) for t in tokens]

    df['psi'] = df['IncLevel1'].apply(
        lambda s: np.mean(vals) if (vals := parse_float_list(s)) else np.nan
    )

    count_cols = [c for c in df.columns if c.startswith(('IJC_SAMPLE', 'SJC_SAMPLE'))]
    def parse_int_sum(values):
        if pd.isna(values): return 0
        total = 0
        for token in str(values).split(','):
            try: total += int(float(token))
            except ValueError: continue
        return total

    for col in count_cols:
        df[col] = df[col].apply(parse_int_sum)
    df['coverage'] = df[count_cols].sum(axis=1)

    return df[(df['coverage'] >= min_coverage) & (~df['psi'].isna())]

# --- Loaders ---

def load_rmats_se(rmats_file, min_coverage=20):
    df = _load_rmats_common(rmats_file, min_coverage)
    if df.empty: return pr.PyRanges()
    
    df = df.rename(columns={'chr': 'Chromosome', 'exonStart_0base': 'Start', 'exonEnd': 'End', 'strand': 'Strand'})
    df[['Start', 'End']] = df[['Start', 'End']].astype(int)
    
    df = _aggregate_exact_duplicates(df)
    return _normalize_strand(df)

def load_rmats_ri(rmats_file, min_coverage=20):
    df = _load_rmats_common(rmats_file, min_coverage)
    if df.empty: return pr.PyRanges()
    
    df['Start'] = df['upstreamEE']
    df['End'] = df['downstreamES']
    df = df.rename(columns={'chr': 'Chromosome', 'strand': 'Strand'})
    df[['Start', 'End']] = df[['Start', 'End']].astype(int)
    
    df = df[df['Start'] < df['End']]
    df = _aggregate_exact_duplicates(df)
    return _normalize_strand(df)

def load_rmats_alt_site(rmats_file, min_coverage=20):
    df = _load_rmats_common(rmats_file, min_coverage)
    if df.empty: return pr.PyRanges()
    
    long_cols = df.rename(columns={'chr': 'Chromosome', 'longExonStart_0base': 'Start', 'longExonEnd': 'End', 'strand': 'Strand'})
    short_cols = df.rename(columns={'chr': 'Chromosome', 'shortES': 'Start', 'shortEE': 'End', 'strand': 'Strand'})
    
    pr_long = _normalize_strand(long_cols[['Chromosome','Start','End','Strand','psi','coverage']])
    pr_short = _normalize_strand(short_cols[['Chromosome','Start','End','Strand','psi','coverage']])
    
    extension = pr_long.subtract(pr_short)
    if len(extension) == 0: return pr.PyRanges()
        
    extension_df = _aggregate_exact_duplicates(extension.df)
    return _normalize_strand(extension_df)

def load_stringtie_gtf(gtf_file):
    return pr.read_gtf(gtf_file)

# --- Generators ---

def generate_hard_labels(stringtie_gr, rmats_se, rmats_ri=None, rmats_a3ss=None, rmats_a5ss=None, psi_low=0.2, psi_high=0.8, min_fragment_len=5):
    # (Matches previous working version)
    stringtie_gr = _ensure_pr(stringtie_gr)
    rmats_se = _ensure_pr(rmats_se)
    rmats_ri = _ensure_pr(rmats_ri)
    rmats_a3ss = _ensure_pr(rmats_a3ss)
    rmats_a5ss = _ensure_pr(rmats_a5ss)

    st_exons = stringtie_gr[stringtie_gr.Feature == 'exon']
    st_transcripts = stringtie_gr[stringtie_gr.Feature == 'transcript']
    
    if len(st_transcripts) == 0:
        st_transcripts = st_exons.merge(by="transcript_id")

    false_exons = rmats_se[rmats_se.psi < psi_low]
    corrected_exons = st_exons.subtract(false_exons) if len(false_exons) else st_exons

    missed_exons = rmats_se[rmats_se.psi > psi_high]
    if len(missed_exons):
        corrected_exons = pr.concat([corrected_exons, missed_exons]).merge()

    if len(rmats_ri):
        ri_high = rmats_ri[rmats_ri.psi > psi_high]
        if len(ri_high):
            ri_high = _normalize_strand(ri_high)
            corrected_exons = pr.concat([corrected_exons, ri_high]).merge()

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

    if min_fragment_len > 0:
        corrected_exons = corrected_exons[corrected_exons.End - corrected_exons.Start >= min_fragment_len]

    introns = st_transcripts.subtract(corrected_exons)

    corrected_exons = _normalize_strand(corrected_exons)
    introns = _normalize_strand(introns)

    corrected_exons.label = 2
    introns.label = 1
    
    return pr.concat([corrected_exons, introns]).sort().df

def generate_soft_labels(stringtie_gr, rmats_se, rmats_ri=None, rmats_a3ss=None, rmats_a5ss=None):
    stringtie_gr = _ensure_pr(stringtie_gr)
    rmats_se = _ensure_pr(rmats_se)
    rmats_ri = _ensure_pr(rmats_ri)
    rmats_a3ss = _ensure_pr(rmats_a3ss)
    rmats_a5ss = _ensure_pr(rmats_a5ss)

    # 1. PREPARE LAYERS
    # -----------------
    
    # Layer: Variables (rMATS)
    vars_list = [gr for gr in (rmats_se, rmats_ri, rmats_a3ss, rmats_a5ss) if len(gr) > 0]
    if vars_list:
        vars_pr = pr.concat(vars_list)
        flat_vars = _resolve_overlaps_mean_psi(vars_pr)
        flat_vars = _normalize_strand(flat_vars)
        if 'psi' not in flat_vars.columns: flat_vars.psi = np.nan
    else:
        flat_vars = _normalize_strand(pr.PyRanges())

    # Layer: Skeleton Exons (PSI=1)
    # Use merge() to get a single flat mask of all exonic regions
    st_exons_raw = stringtie_gr[stringtie_gr.Feature == 'exon']
    st_exons_flat = _normalize_strand(st_exons_raw.merge())
    st_exons_flat.psi = 1.0

    # Layer: Canvas (Transcripts)
    st_transcripts_raw = stringtie_gr[stringtie_gr.Feature == 'transcript']
    if len(st_transcripts_raw) == 0:
        st_transcripts_raw = stringtie_gr[stringtie_gr.Feature == 'exon'].merge(by="transcript_id")
    # This canvas represents the universe of all transcribed regions
    canvas = _normalize_strand(st_transcripts_raw.merge())

    # 2. STAMPING (Strict Priority Order)
    # -----------------------------------
    final_parts = []
    
    # Priority A: Variable Regions (rMATS)
    if len(flat_vars) > 0:
        # Keep only vars that fall within the gene models
        on_canvas = flat_vars.intersect(canvas)
        
        if len(on_canvas) > 0:
            on_canvas = _normalize_strand(on_canvas)
            final_parts.append(on_canvas)
            
            # CRITICAL: Remove these regions from the canvas
            # We must aggressively normalize after every subtract/intersect
            canvas = canvas.subtract(on_canvas)
            canvas = _normalize_strand(canvas)

    # Priority B: Constitutive Exons (PSI=1.0)
    # We stamp PSI=1 only if the region is in st_exons_flat AND wasn't already covered by rMATS
    if len(st_exons_flat) > 0 and len(canvas) > 0:
        # Intersect the REMAINING canvas with the exon mask
        # This gives us parts of the canvas that are exonic but not variable
        const_exons = canvas.intersect(st_exons_flat)
        
        if len(const_exons) > 0:
            const_exons = _normalize_strand(const_exons)
            const_exons.psi = 1.0
            final_parts.append(const_exons)
            
            # Remove from canvas
            canvas = canvas.subtract(const_exons)
            canvas = _normalize_strand(canvas)

    # Priority C: Introns (PSI=0.0)
    # Anything left on the canvas is, by definition, an Intron
    if len(canvas) > 0:
        canvas.psi = 0.0
        final_parts.append(canvas)

    # 3. ASSEMBLY
    # -----------
    if not final_parts:
        return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'p_exon', 'Strand'])

    # Explicitly ensure strand column exists on all parts before concat
    # (Double check to prevent that previous AttributeError)
    safe_parts = []
    for p in final_parts:
        if 'Strand' not in p.df.columns:
            p = _normalize_strand(p)
        safe_parts.append(p)

    combined = pr.concat(safe_parts).sort()
    df_out = _aggregate_exact_duplicates(combined.df)
    
    return df_out[['Chromosome', 'Start', 'End', 'psi', 'Strand']].rename(columns={'psi': 'p_exon'})