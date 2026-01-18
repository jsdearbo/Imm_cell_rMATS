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

def _load_rmats_common(rmats_file, min_coverage):
    try:
        df = pd.read_csv(rmats_file, sep='\t')
    except FileNotFoundError:
        print(f"[WARN] File not found: {rmats_file}")
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        print(f"[WARN] File empty: {rmats_file}")
        return pd.DataFrame()
    except pd.errors.ParserError as exc:
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
    
    # We must construct clean PyRanges for subtraction
    long_cols = df.rename(columns={'chr': 'Chromosome', 'longExonStart_0base': 'Start', 'longExonEnd': 'End', 'strand': 'Strand'})
    short_cols = df.rename(columns={'chr': 'Chromosome', 'shortES': 'Start', 'shortEE': 'End', 'strand': 'Strand'})
    
    # Create temp PyRanges for the operation
    pr_long = _normalize_strand(long_cols[['Chromosome','Start','End','Strand','psi','coverage']])
    pr_short = _normalize_strand(short_cols[['Chromosome','Start','End','Strand','psi','coverage']])
    
    # Calculate Extension
    extension = pr_long.subtract(pr_short)
    
    if len(extension) == 0:
        return pr.PyRanges()
        
    extension_df = _aggregate_exact_duplicates(extension.df)
    return _normalize_strand(extension_df)

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

    # 1. Collect Variable Regions
    vars_list = [gr for gr in (rmats_se, rmats_ri, rmats_a3ss, rmats_a5ss) if len(gr) > 0]
    
    vars_pr = pr.concat(vars_list) if vars_list else pr.PyRanges()
    
    # If we have variables, we normalize and de-duplicate them immediately
    if len(vars_pr):
        vars_pr = _normalize_strand(pr.PyRanges(_aggregate_exact_duplicates(vars_pr.df)))

    # 2. Constitutive Exons
    st_exons = stringtie_gr[stringtie_gr.Feature == 'exon']
    const_exons = st_exons.subtract(vars_pr) if len(vars_pr) else st_exons
    const_exons.psi = 1.0
    
    # 3. Constitutive Introns
    st_transcripts = stringtie_gr[stringtie_gr.Feature == 'transcript']
    if len(st_transcripts) == 0:
        st_transcripts = st_exons.merge(by="transcript_id")

    known_features = pr.concat([const_exons, vars_pr]) if len(vars_pr) else const_exons
    const_introns = st_transcripts.subtract(known_features)
    const_introns.psi = 0.0
    
    # 4. Combine
    combined_parts = [const_exons, const_introns]
    if len(vars_pr):
        combined_parts.append(vars_pr)
        
    combined = pr.concat(combined_parts).sort()
    
    # Final cleanup of any duplicate boundaries
    df_out = _aggregate_exact_duplicates(combined.df)
    
    return df_out[['Chromosome', 'Start', 'End', 'psi']].rename(columns={'psi': 'p_exon'})