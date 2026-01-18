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
    Resolves EXACT coordinate duplicates with sensible defaults.
    Prefers rows carrying an explicit priority flag (when present), otherwise
    aggregates numeric values such as psi/coverage.
    """
    if df.empty:
        return df

    df = df.reset_index(drop=True)

    group_cols = ['Chromosome', 'Start', 'End']
    if 'Strand' in df.columns:
        group_cols.append('Strand')

    priority_col = None
    for candidate in ('source_priority', 'priority'):
        if candidate in df.columns:
            priority_col = candidate
            break

    if priority_col is not None:
        sort_cols = group_cols + [priority_col]
        df = df.sort_values(sort_cols)
        dedup = df.drop_duplicates(subset=group_cols, keep='first')
        return dedup.drop(columns=[c for c in (priority_col, 'source') if c in dedup.columns])

    agg_rules = {}
    if 'psi' in df.columns:
        agg_rules['psi'] = 'mean'
    if 'coverage' in df.columns:
        agg_rules['coverage'] = 'sum'

    if agg_rules:
        return df.groupby(group_cols, observed=True, as_index=False).agg(agg_rules)

    return df.groupby(group_cols, observed=True, as_index=False).first()

def _resolve_overlaps_mean_psi(gr):
    """Split overlaps into atomic intervals and average PSI over overlaps."""
    if gr is None or len(gr) == 0:
        return pr.PyRanges()

    df = gr.df.reset_index(drop=True)
    results = []

    for (chrom, strand), group in df.groupby(['Chromosome', 'Strand'], observed=True):
        if group.empty:
            continue

        boundaries = np.unique(np.concatenate([group.Start.values, group.End.values]))
        if len(boundaries) < 2:
            continue

        starts = boundaries[:-1]
        ends = boundaries[1:]
        valid = starts < ends
        starts = starts[valid]
        ends = ends[valid]
        if len(starts) == 0:
            continue

        elem_df = pd.DataFrame({
            'Chromosome': chrom,
            'Start': starts,
            'End': ends,
            'Strand': strand
        })
        elem_pr = pr.PyRanges(elem_df)
        group_pr = pr.PyRanges(group)
        joined = elem_pr.join(group_pr)
        if joined is None or len(joined) == 0:
            continue

        agg = (
            joined.df
            .groupby(['Chromosome', 'Start', 'End', 'Strand'], observed=True)['psi']
            .mean()
            .reset_index()
        )
        results.append(agg)

    if not results:
        return pr.PyRanges()

    final_df = pd.concat(results, ignore_index=True)
    return _normalize_strand(pr.PyRanges(final_df))

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

    # Variable regions (rMATS) flattened to disjoint PSI-bearing segments
    vars_list = [gr for gr in (rmats_se, rmats_ri, rmats_a3ss, rmats_a5ss) if len(gr) > 0]
    if vars_list:
        vars_pr = pr.concat(vars_list)
        flat_vars = _resolve_overlaps_mean_psi(vars_pr)
        flat_vars = _normalize_strand(flat_vars) if len(flat_vars) else pr.PyRanges()
        if len(flat_vars) and 'psi' not in flat_vars.columns:
            flat_vars.psi = np.nan
    else:
        flat_vars = pr.PyRanges()

    # Constitutive exons flattened from StringTie
    st_exons = stringtie_gr[stringtie_gr.Feature == 'exon']
    flat_exons = st_exons.merge()
    flat_exons = _normalize_strand(flat_exons) if len(flat_exons) else pr.PyRanges()

    # Transcript canvases ensure continuous coverage across loci
    st_transcripts = stringtie_gr[stringtie_gr.Feature == 'transcript']
    if len(st_transcripts) == 0:
        st_transcripts = st_exons.merge(by="transcript_id")
    flat_transcripts = st_transcripts.merge()
    flat_transcripts = _normalize_strand(flat_transcripts) if len(flat_transcripts) else pr.PyRanges()

    if len(flat_transcripts) == 0:
        return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'p_exon', 'Strand'])

    vars_df = flat_vars.df if len(flat_vars) else pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Strand', 'psi'])
    exons_df = flat_exons.df if len(flat_exons) else pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Strand'])
    tx_df = flat_transcripts.df

    results = []
    priority_map = {'variable': 0, 'exon': 1, 'intron': 2}

    for _, tx in tx_df.iterrows():
        chrom = tx['Chromosome']
        tx_start = int(tx['Start'])
        tx_end = int(tx['End'])
        strand = tx.get('Strand', '.')

        boundaries = {tx_start, tx_end}

        exon_slice = pd.DataFrame(columns=['Start', 'End'])
        if not exons_df.empty:
            mask = (
                (exons_df['Chromosome'] == chrom) &
                (exons_df['Strand'] == strand) &
                (exons_df['End'] > tx_start) &
                (exons_df['Start'] < tx_end)
            )
            exon_slice = exons_df[mask].copy()
            if not exon_slice.empty:
                exon_slice['Start'] = exon_slice['Start'].clip(lower=tx_start)
                exon_slice['End'] = exon_slice['End'].clip(upper=tx_end)
                boundaries.update(exon_slice['Start'].astype(int).tolist())
                boundaries.update(exon_slice['End'].astype(int).tolist())

        var_slice = pd.DataFrame(columns=['Start', 'End', 'psi'])
        if not vars_df.empty:
            mask = (
                (vars_df['Chromosome'] == chrom) &
                (vars_df['Strand'] == strand) &
                (vars_df['End'] > tx_start) &
                (vars_df['Start'] < tx_end)
            )
            var_slice = vars_df[mask].copy()
            if not var_slice.empty:
                var_slice['Start'] = var_slice['Start'].clip(lower=tx_start)
                var_slice['End'] = var_slice['End'].clip(upper=tx_end)
                boundaries.update(var_slice['Start'].astype(int).tolist())
                boundaries.update(var_slice['End'].astype(int).tolist())

        if len(boundaries) <= 1:
            continue

        cuts = np.array(sorted(boundaries), dtype=int)
        seg_starts = cuts[:-1]
        seg_ends = cuts[1:]

        for seg_start, seg_end in zip(seg_starts, seg_ends):
            if seg_end <= seg_start:
                continue

            psi_val = None
            source = 'intron'

            if not var_slice.empty:
                var_hit = var_slice[
                    (var_slice['Start'] <= seg_start) &
                    (var_slice['End'] >= seg_end)
                ]
                if not var_hit.empty:
                    psi_candidate = var_hit.iloc[0]['psi']
                    if not pd.isna(psi_candidate):
                        psi_val = float(psi_candidate)
                        source = 'variable'

            if psi_val is None and not exon_slice.empty:
                exon_hit = exon_slice[
                    (exon_slice['Start'] <= seg_start) &
                    (exon_slice['End'] >= seg_end)
                ]
                if not exon_hit.empty:
                    psi_val = 1.0
                    source = 'exon'

            if psi_val is None:
                psi_val = 0.0
                source = 'intron'

            results.append({
                'Chromosome': chrom,
                'Start': int(seg_start),
                'End': int(seg_end),
                'Strand': strand if strand in ['+', '-'] else '.',
                'psi': psi_val,
                'source_priority': priority_map[source]
            })

    if not results:
        return pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'p_exon', 'Strand'])

    result_df = pd.DataFrame(results)
    sort_cols = ['Chromosome', 'Start', 'End', 'Strand', 'source_priority']
    existing_sort_cols = [c for c in sort_cols if c in result_df.columns]
    result_df = result_df.sort_values(existing_sort_cols)

    group_cols = ['Chromosome', 'Start', 'End', 'Strand']
    existing_group_cols = [c for c in group_cols if c in result_df.columns]
    if existing_group_cols:
        result_df = (
            result_df
            .groupby(existing_group_cols, as_index=False)
            .first()
        )

    result_df = result_df.drop(columns=['source_priority'], errors='ignore')

    result_pr = _normalize_strand(result_df)
    result_df = result_pr.sort().df

    if 'Strand' not in result_df.columns:
        result_df['Strand'] = '.'

    cols = ['Chromosome', 'Start', 'End', 'psi', 'Strand']
    return result_df[cols].rename(columns={'psi': 'p_exon'})