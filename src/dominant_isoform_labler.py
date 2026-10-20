import pandas as pd
import gtfparse
import os
import glob

# =====================
# Configuration
# =====================
GTF_DIR = "/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_per_rep"
MERGED_DIR = f"{GTF_DIR}/merged"
REQUANT_DIR = f"{GTF_DIR}/requantified"
DOMINANT_GTF_DIR = f"{GTF_DIR}/dominant_isoforms"

MIN_TRANSCRIPT_TPM = 1.0
MIN_TRANSCRIPT_COV = 2.0
MIN_EXON_LEN = 0

CELL_TYPES = ['B_Cell', 'T_Cell', 'Macrophage']

os.makedirs(DOMINANT_GTF_DIR, exist_ok=True)

print("[INFO] Building dominant transcript GTFs (no binning)")
print(f"[INFO] Filters: mean_TPM>={MIN_TRANSCRIPT_TPM}, mean_cov>={MIN_TRANSCRIPT_COV}, exon_len>={MIN_EXON_LEN}")

# =====================
# Helpers
# =====================
def _safe_to_numeric(df, col):
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    else:
        df[col] = pd.NA


def _collect_replicate_gtfs(cell_type):
    """Find all re-quantified GTFs for a cell type."""
    pattern = os.path.join(REQUANT_DIR, f"{cell_type}_*", f"{cell_type}_*.gtf")
    paths = sorted(glob.glob(pattern))
    return paths


def get_dominant_tx_ids(cell_type):
    """
    Read re-quantified per-replicate GTFs, average TPM/cov across replicates,
    then select the dominant (highest mean-TPM) transcript per reference gene.
    """
    print(f"[INFO] Processing {cell_type}")

    rep_gtfs = _collect_replicate_gtfs(cell_type)
    if not rep_gtfs:
        print(f"[WARN] No re-quantified GTFs found for {cell_type}")
        return [], None
    print(f"\t- Found {len(rep_gtfs)} replicate GTFs")

    # ---- Read transcript-level TPM/cov from each replicate ----
    rep_frames = []
    for i, gpath in enumerate(rep_gtfs):
        sample = os.path.basename(os.path.dirname(gpath))
        df = gtfparse.read_gtf(gpath).to_pandas()
        _safe_to_numeric(df, 'TPM')
        _safe_to_numeric(df, 'cov')
        tx = df[df['feature'] == 'transcript'][['transcript_id', 'gene_id', 'TPM', 'cov']].copy()
        # Carry ref_gene_id if present
        if 'ref_gene_id' in df.columns:
            tx['ref_gene_id'] = df.loc[tx.index, 'ref_gene_id']
        tx = tx.rename(columns={'TPM': f'TPM_{i}', 'cov': f'cov_{i}'})
        rep_frames.append(tx)

    # ---- Merge replicates on transcript_id and compute means ----
    merged = rep_frames[0]
    for frame in rep_frames[1:]:
        # Keep only quant columns for the join
        cols_to_join = ['transcript_id'] + [c for c in frame.columns if c.startswith(('TPM_', 'cov_'))]
        merged = merged.merge(frame[cols_to_join], on='transcript_id', how='outer')

    tpm_cols = [c for c in merged.columns if c.startswith('TPM_')]
    cov_cols = [c for c in merged.columns if c.startswith('cov_')]
    merged['mean_TPM'] = merged[tpm_cols].mean(axis=1)
    merged['mean_cov'] = merged[cov_cols].mean(axis=1)

    n_total = len(merged)
    print(f"\t- Total transcripts across replicates: {n_total}")

    # ---- Filter ----
    merged = merged[merged['mean_TPM'].fillna(0) >= MIN_TRANSCRIPT_TPM]
    merged = merged[merged['mean_cov'].fillna(0) >= MIN_TRANSCRIPT_COV]

    if merged.empty:
        print(f"[WARN] No transcripts passed filters for {cell_type}")
        return [], None

    print(f"\t- Transcripts passing TPM/cov filters: {len(merged)}")

    # ---- Group by reference gene ----
    if 'ref_gene_id' in merged.columns and merged['ref_gene_id'].notna().any():
        group_col = 'ref_gene_id'
    elif 'gene_id' in merged.columns:
        group_col = 'gene_id'
    else:
        raise ValueError(f"No gene grouping column found for {cell_type}")
    print(f"\t- Grouping by: {group_col}")

    dominant = merged.loc[merged.groupby(group_col, observed=True)['mean_TPM'].idxmax()]
    dominant_tx_ids = dominant['transcript_id'].tolist()
    print(f"\t- Dominant transcripts before exon-length filter: {len(dominant_tx_ids)}")

    # ---- Exon-length filter (read from first replicate for structure) ----
    struct_df = gtfparse.read_gtf(rep_gtfs[0]).to_pandas()
    exon_df = struct_df[
        (struct_df['feature'] == 'exon') &
        (struct_df['transcript_id'].isin(dominant_tx_ids))
    ][['transcript_id', 'start', 'end']].copy()

    if exon_df.empty:
        print(f"[WARN] No exon entries found for dominant transcripts in {cell_type}")
        return [], None

    exon_df['exon_len'] = exon_df['end'] - exon_df['start'] + 1
    valid_tx = set(exon_df[exon_df['exon_len'] >= MIN_EXON_LEN]['transcript_id'])
    filtered_tx_ids = [tx for tx in dominant_tx_ids if tx in valid_tx]

    print(f"\t- Dominant transcripts after exon-length filter: {len(filtered_tx_ids)}")

    # Return the merged-model GTF path for writing (has clean structure)
    merged_gtf = os.path.join(MERGED_DIR, f"{cell_type}_merged.gtf")
    # Fall back to the first replicate GTF if merged model doesn't exist
    src_gtf = merged_gtf if os.path.exists(merged_gtf) else rep_gtfs[0]
    return filtered_tx_ids, src_gtf


def write_dominant_gtf(src_path, keep_tx_ids, dest_path, min_exon_len):
    keep = set(keep_tx_ids)
    written = 0
    with open(src_path) as inp, open(dest_path, "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if 'transcript_id "' not in line:
                continue
            tx = line.split('transcript_id "')[1].split('"')[0]
            if tx not in keep:
                continue

            feature = fields[2]
            if feature == 'exon' and min_exon_len > 0:
                start = int(fields[3])
                end = int(fields[4])
                if (end - start + 1) < min_exon_len:
                    continue

            if feature in {'transcript', 'exon', 'five_prime_utr', 'three_prime_utr', 'CDS', 'start_codon', 'stop_codon'}:
                out.write(line)
                written += 1
    return written


# =====================
# Build dominant GTF per cell type
# =====================

for ct in CELL_TYPES:
    tx_ids, src_gtf = get_dominant_tx_ids(ct)
    if not tx_ids:
        print(f"[WARN] No dominant transcripts kept for {ct}, skipping write.")
        continue

    dom_gtf = os.path.join(DOMINANT_GTF_DIR, f"{ct}_dominant.gtf")
    n_lines = write_dominant_gtf(src_gtf, tx_ids, dom_gtf, MIN_EXON_LEN)
    print(f"[INFO] Wrote {dom_gtf} ({n_lines} lines from {os.path.basename(src_gtf)})")

print(f"[DONE] Dominant GTFs written to {DOMINANT_GTF_DIR}")
