import pandas as pd
import numpy as np
import gtfparse
import os

# Set configuration parameters
BIN_SIZE = 32 # The fixed genomic resolution for the model
GTF_DIR = "/mnt/netfiles/jsdearbo/scratch/Imm_cell_rMATS/results/stringtie_out"
OUTPUT_DIR = f"{GTF_DIR}/multi_task_targets"
DOMINANT_GTF_DIR = f"{GTF_DIR}/dominant_isoforms"
REFERENCE_GTF = "/mnt/netfiles/genomes/STAR/mm10/gencode.vM23.annotation.gtf" # Used for defining the genome coordinates

# Define the list of cell types and their aggregated GTF files
CELL_TYPES = ['B_Cell', 'T_Cell', 'Macrophage']
GTF_PATHS = { 
    ct: os.path.join(GTF_DIR, f'{ct}_aggregate.gtf') 
    for ct in CELL_TYPES
}

# Class mapping
CLASS_MAP = {
    'exon': 2, 
    'UTR': 3, 
    'intron': 1, 
    'desert': 0 
}

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(DOMINANT_GTF_DIR, exist_ok=True)
print(f"Target bins size: {BIN_SIZE} bp")

def generate_bins(gtf_path, bin_size):
    """Reads GTF to infer max chrom length and generate fixed genomic bins."""
    
    print("\n[INFO] Loading reference GTF to define genome boundaries...")
    # Load only necessary columns for efficiency
    ref_df = gtfparse.read_gtf(gtf_path).to_pandas()[['seqname', 'end']]
    
    # Find the maximum end coordinate for each chromosome
    chrom_lengths = ref_df.groupby('seqname')['end'].max()
    
    all_bins = []
    for chrom, max_len in chrom_lengths.items():
        start_coords = np.arange(1, max_len + 1, bin_size)
        end_coords = np.minimum(start_coords + bin_size - 1, max_len)
        
        chrom_bins = pd.DataFrame({
            'chrom': chrom,
            'start': start_coords,
            'end': end_coords,
        })
        all_bins.append(chrom_bins)
        
    master_bins_df = pd.concat(all_bins).reset_index(drop=True)
    print(f"Generated {len(master_bins_df)} fixed bins of size {bin_size} bp.")
    return master_bins_df

# Run bin generation
master_bins_df = generate_bins(REFERENCE_GTF, BIN_SIZE)
master_bins_df['bin_id'] = master_bins_df.index
master_bins_df.head()

import pyranges as pr
def get_dominant_features(gtf_path, cell_type):
    """Loads GTF, finds dominant transcript (highest TPM) per gene, and derives features."""
    print(f"[INFO] Processing {cell_type} GTF: {gtf_path}")
    
    # Load all records (gtfparse handles the quotes and types)
    full_gtf = gtfparse.read_gtf(gtf_path).to_pandas()
    full_gtf['TPM'] = pd.to_numeric(full_gtf['TPM'], errors='coerce')
    
    # Filter for transcripts to identify the dominant ones (using the 'transcript' feature)
    tx_df = full_gtf[
        (full_gtf['feature'] == 'transcript') &
        (full_gtf['TPM'] > 0) 
    ].copy()
    
    # Ensure TPM is numeric for sorting
    tx_df['TPM'] = pd.to_numeric(tx_df['TPM'], errors='coerce')
    
    # 1. Find the dominant transcript (highest TPM) per gene
    dominant_tx = tx_df.loc[tx_df.groupby('gene_id')['TPM'].idxmax()]
    print(f"\t- Identified {len(dominant_tx)} dominant transcripts.")
    
    # List of dominant transcript IDs
    dominant_tx_ids = dominant_tx['transcript_id'].tolist()
    
    # 2. Extract features (exons and UTRs) belonging to the dominant transcripts
    dominant_features_df = full_gtf[
        full_gtf['transcript_id'].isin(dominant_tx_ids) &
        full_gtf['feature'].isin(['exon', 'five_prime_utr', 'three_prime_utr'])
    ].copy()
    
    # 3. Simplify feature type for labeling
    dominant_features_df['class_feature'] = dominant_features_df['feature'].map({
        'exon': 'exon', 
        'five_prime_utr': 'UTR', 
        'three_prime_utr': 'UTR' 
    })
    
    # 4. Derive Intron features based on the dominant transcripts
    
    # Sort features by transcript and start position (crucial for deriving gaps)
    sorted_exons = dominant_features_df[
        dominant_features_df['feature'] == 'exon'
    ].sort_values(by=['transcript_id', 'start'])
    
    intron_list = []
    for tx_id, group in sorted_exons.groupby('transcript_id'):
        # Get the start of the next exon (shifted up)
        next_exon_start = group['start'].shift(-1)
        # Get the end of the current exon
        current_exon_end = group['end']
        
        # Intron start = current_exon_end + 1
        # Intron end = next_exon_start - 1
        intron_starts = current_exon_end + 1
        intron_ends = next_exon_start - 1
        
        # Create intron records where the gap is valid (start < end)
        intron_records = group.loc[intron_starts < intron_ends].copy()
        intron_records['start'] = intron_starts
        intron_records['end'] = intron_ends
        intron_records['class_feature'] = 'intron'
        intron_list.append(intron_records)

    if intron_list:
        intron_df = pd.concat(intron_list)
        dominant_features_df = pd.concat([dominant_features_df, intron_df], ignore_index=True)
        print(f"\t- Derived {len(intron_df)} intron regions.")
    
    # Final cleaning and selecting columns needed for binning
    dominant_features_df = dominant_features_df[[
        'seqname', 'start', 'end', 'class_feature', 'transcript_id'
    ]].rename(columns={'seqname': 'chrom'}).drop_duplicates()
    
    print(f"\t- Total unique features to map: {len(dominant_features_df)}")
    return dominant_features_df, dominant_tx_ids

def write_dominant_gtf(src_path, transcript_ids, dest_path):
    tag = 'transcript_id "'
    keep = set(transcript_ids)
    with open(src_path, "r") as inp, open(dest_path, "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            idx = line.find(tag)
            if idx == -1:
                continue
            start = idx + len(tag)
            end = line.find('"', start)
            tx_id = line[start:end]
            if tx_id in keep:
                out.write(line)

master_labels = master_bins_df[['bin_id', 'chrom', 'start', 'end']].copy()

# Global priority (used to resolve overlaps, higher value wins)
PRIORITY_MAP = {'exon': 3, 'UTR': 2, 'intron': 1, 'desert': 0}

for cell_type in CELL_TYPES:
    gtf_path = GTF_PATHS[cell_type]
    
    if not os.path.exists(gtf_path):
        print(f"[ERROR] File not found: {gtf_path}. Skipping {cell_type}.")
        continue
        
    # Get the features for this cell type
    features_df, dominant_tx_ids = get_dominant_features(gtf_path, cell_type)
    
    # Write the dominant isoform GTF for this cell type
    dominant_gtf_path = os.path.join(DOMINANT_GTF_DIR, f"{cell_type}_dominant.gtf")
    write_dominant_gtf(gtf_path, dominant_tx_ids, dominant_gtf_path)
    
    # --- Binning and Labeling --- 
    
    # Assign priority for overlap resolution
    features_df['priority'] = features_df['class_feature'].map(PRIORITY_MAP)
    features_df['class'] = features_df['class_feature'].map(CLASS_MAP)
    
    # Perform a genomic merge (overlap check)
    # This relies on the core idea: for each bin, find all overlapping features.
    
    # Merge features with master bins based on coordinate overlap
    # 1. Match chromosome
    # 2. Bin_Start <= Feature_End
    # 3. Bin_End >= Feature_Start
    bin_gr = pr.PyRanges(master_bins_df.rename(columns={"start": "Start", "end": "End", "chrom": "Chromosome"}))
    clean_features = features_df.dropna(subset=["start", "end"]).copy()
    clean_features[["start", "end"]] = clean_features[["start", "end"]].astype(int)

    feat_gr = pr.PyRanges(clean_features.rename(
        columns={"start": "Start", "end": "End", "chrom": "Chromosome"}
    ))

    overlaps = bin_gr.join(feat_gr, how="left").as_df()
    overlaps.rename(columns={"Start": "start_bin", "End": "end_bin"}, inplace=True)
    overlaps = overlaps.dropna(subset=["priority"])  # drop desert bins
    
    # Handle Overlap Resolution (Prioritization)
    # Group by bin_id and keep the feature with the highest priority
    max_priority = overlaps.loc[overlaps.groupby('bin_id')['priority'].idxmax()]
    
    # Reduce to just the bin_id and the assigned class
    final_labels = max_priority[['bin_id', 'class']].drop_duplicates()
    
    # Rename the class column for the multi-task matrix
    label_col_name = f'{cell_type}_Class'
    final_labels = final_labels.rename(columns={'class': label_col_name})
    
    # Merge back into the master label table (using left merge)
    master_labels = master_labels.merge(final_labels, on='bin_id', how='left')
    
    # Fill NaN values (where no feature overlapped, or Desert) with Class 0
    master_labels[label_col_name] = master_labels[label_col_name].fillna(CLASS_MAP['desert']).astype(int)
    
    print(f"[INFO] {cell_type} labeling complete.")
    print("------------------")
    

# Final cleanup and export
master_labels = master_labels.drop(columns='bin_id')
output_path = os.path.join(OUTPUT_DIR, "multi_task_classification_targets.tsv")
master_labels.to_csv(output_path, sep='\t', index=False)

print(f"\nFINAL TARGET MATRIX SAVED to {output_path}.")
print(f"Sample of the Multi-Task Classification Matrix:\n{master_labels.head()}")