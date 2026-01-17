import pandas as pd
import numpy as np
import pyranges as pr

def load_rmats_se(rmats_file, min_coverage=20):
    """
    Loads rMATS Skipped Exon (SE) file.
    Returns PyRanges with 'psi' column.
    """
    df = pd.read_csv(rmats_file, sep='\t')
    
    # Helper to parse PSI lists (e.g., "0.1,0.2,0.15")
    def parse_float_list(values):
        tokens = [t for t in str(values).split(',') if t not in {'', 'NA'}]
        return [float(t) for t in tokens]

    # Calculate mean PSI from IncLevel1 (assuming replicate 1 is the target cell type)
    df['psi'] = df['IncLevel1'].apply(
        lambda s: np.mean(vals) if (vals := parse_float_list(s)) else np.nan
    )

    # Calculate coverage
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

    # Filter
    df = df[(df['coverage'] >= min_coverage) & (~df['psi'].isna())]

    # Convert to standard PyRanges format
    # rMATS SE: exonStart_0base, exonEnd
    pr_df = df[['chr', 'exonStart_0base', 'exonEnd', 'strand', 'psi']].copy()
    pr_df = pr_df.rename(columns={
        'chr': 'Chromosome',
        'exonStart_0base': 'Start',
        'exonEnd': 'End',
        'strand': 'Strand'
    })
    # Ensure standard types
    pr_df[['Start', 'End']] = pr_df[['Start', 'End']].astype(int)

    return pr.PyRanges(pr_df)

def load_rmats_ri(rmats_file, min_coverage=20):
    """
    Loads rMATS Retained Intron (RI) file.
    Returns PyRanges of the *INTRON* region with 'psi' column.
    """
    try:
        df = pd.read_csv(rmats_file, sep='\t')
    except FileNotFoundError:
        print(f"[WARN] RI File not found: {rmats_file}. Returning empty PyRanges.")
        return pr.PyRanges()

    # Reuse parsing logic
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
    df = df[(df['coverage'] >= min_coverage) & (~df['psi'].isna())]

    # Crucial: Define the Intron Region.
    # rMATS RI: upstreamEE (end of upstream exon) to downstreamES (start of downstream exon)
    # rMATS coordinates are often 0-based/1-based mixed.
    # upstreamEE is 1-based end of upstream? 
    # Usually: Intron Start = upstreamEE, Intron End = downstreamES (0-based)
    # Let's map strictly to the gap.
    
    df['Start'] = df['upstreamEE']
    df['End'] = df['downstreamES']
    
    pr_df = df[['chr', 'Start', 'End', 'strand', 'psi']].copy()
    pr_df = pr_df.rename(columns={'chr': 'Chromosome', 'strand': 'Strand'})
    pr_df[['Start', 'End']] = pr_df[['Start', 'End']].astype(int)
    
    # Filter out bad coordinates (e.g. if overlap causes negative length)
    pr_df = pr_df[pr_df['Start'] < pr_df['End']]

    return pr.PyRanges(pr_df)

def load_stringtie_gtf(gtf_file):
    """
    Loads StringTie GTF. Returns ALL features (exons and transcripts).
    """
    gr = pr.read_gtf(gtf_file)
    return gr

def generate_hard_labels(stringtie_gr, rmats_se, rmats_ri=None, psi_low=0.2, psi_high=0.8):
    """
    Generates Hard Labels: Intron=1, Exon=2.
    1. Correct StringTie Exons using rMATS SE (Remove Low PSI, Add High PSI).
    2. Add High PSI Retained Introns to Exons.
    3. Define Introns as (Transcript Extent - Exons).
    """
    
    # 1. Separate Features
    st_exons = stringtie_gr[stringtie_gr.Feature == 'exon']
    st_transcripts = stringtie_gr[stringtie_gr.Feature == 'transcript']
    
    # If transcripts missing, infer from exons (group by transcript_id)
    if len(st_transcripts) == 0:
        st_transcripts = st_exons.merge(by="transcript_id")

    # 2. Correct Exons (SE)
    # Remove False Positives (StringTie says Exon, rMATS says PSI < low)
    false_exons = rmats_se[rmats_se.psi < psi_low]
    # Keep ST exons that do NOT overlap false_exons
    corrected_exons = st_exons.subtract(false_exons)

    # Add False Negatives (StringTie missed it, rMATS says PSI > high)
    missed_exons = rmats_se[rmats_se.psi > psi_high]
    if len(missed_exons) > 0:
        # FIX: Use concat to ADD, not join (which intersects)
        missed_exons = _normalize_strand(missed_exons)
        corrected_exons = pr.concat([corrected_exons, missed_exons]).merge()

    # 3. Handle Retained Introns (RI)
    # If RI PSI > high, it behaves like an Exon (Sequence is present)
    if rmats_ri is not None and len(rmats_ri) > 0:
        retained_introns = rmats_ri[rmats_ri.psi > psi_high]
        if len(retained_introns) > 0:
            # FIX: Use concat to ADD
            retained_introns = _normalize_strand(retained_introns)
            corrected_exons = pr.concat([corrected_exons, retained_introns]).merge()

    # 4. Define Introns (The Subtraction Step)
    # Introns = Transcript Body - Corrected Exons
    # This automatically handles the gaps created/filled by SE/RI corrections
    introns = st_transcripts.subtract(corrected_exons)

    # Normalize strands to prevent PyRanges mismatch errors during concatenation
    corrected_exons = _normalize_strand(corrected_exons)
    introns = _normalize_strand(introns)

    # 5. Label Assignment
    corrected_exons.label = 2
    introns.label = 1
    
    # Combine
    combined = pr.concat([corrected_exons, introns]).sort()
    
    return combined.df

def generate_soft_labels(stringtie_gr, rmats_se, rmats_ri=None):
    """
    Generates Soft Labels: Score = P(Exon).
    Exon = 1.0 (or SE PSI)
    Intron = 0.0 (or RI PSI)
    """
    
    st_exons = stringtie_gr[stringtie_gr.Feature == 'exon']
    st_transcripts = stringtie_gr[stringtie_gr.Feature == 'transcript']
    
    if len(st_transcripts) == 0:
        st_transcripts = st_exons.merge(by="transcript_id")

    # 1. Define Regions
    
    # Variable Exons (SE) -> Use PSI
    vars_se = rmats_se[['Chromosome', 'Start', 'End', 'psi']]
    
    # Variable Introns (RI) -> Use PSI (Retained Intron PSI is effectively P(Exon))
    vars_ri = pr.PyRanges()
    if rmats_ri is not None and len(rmats_ri) > 0:
        vars_ri = rmats_ri[['Chromosome', 'Start', 'End', 'psi']]

    # Constitutive Exons: StringTie Exons NOT overlapping SE or RI
    # (Subtract variable regions to avoid double counting)
    const_exons = st_exons.subtract(vars_se).subtract(vars_ri)
    const_exons.psi = 1.0
    
    # 2. Constitutive Introns
    # Everything in transcript that is NOT an Exon or Variable Region
    # P(Exon) = 0.0
    known_features = pr.concat([const_exons, vars_se, vars_ri])
    const_introns = st_transcripts.subtract(known_features)
    const_introns.psi = 0.0
    
    # 3. Combine All
    # Order matters less for concat, but let's be clean
    combined = pr.concat([const_exons, const_introns, vars_se, vars_ri]).sort()
    
    df_out = combined.df[['Chromosome', 'Start', 'End', 'psi']].copy()
    df_out = df_out.rename(columns={'psi': 'p_exon'})
    
    return df_out

def _normalize_strand(gr):
    df = gr.df.copy()
    if "Strand" not in df.columns:
        df["Strand"] = "."
    df["Strand"] = df["Strand"].astype("category")
    return pr.PyRanges(df)