import pandas as pd
import numpy as np
import pyranges as pr
import os
import argparse
import gtfparse

# convert bed files with labels into dense bin labels
# for both hard and soft labeling schemes

# --- 1. Reusing your Binning Logic ---
def generate_bins(gtf_path, bin_size, shift=0):
    """Generates the master bin dataframe."""
    print(f"[INFO] Loading reference GTF: {gtf_path}")
    # Using a simplified reader to avoid gtfparse version dependency issues if unrelated
    # Assuming standard GTF/GFF format
    ref_df = gtfparse.read_gtf(gtf_path).to_pandas()[['seqname', 'end']]
    chrom_lengths = ref_df.groupby('seqname')['end'].max()
    
    all_bins = []
    for chrom, max_len in chrom_lengths.items():
        start_coords = np.arange(1 + shift, max_len + 1, bin_size)
        end_coords = np.minimum(start_coords + bin_size - 1, max_len)
        chrom_bins = pd.DataFrame({
            'chrom': chrom,
            'start': start_coords,
            'end': end_coords,
        })
        # Filter valid
        chrom_bins = chrom_bins[chrom_bins['start'] < chrom_bins['end']]
        all_bins.append(chrom_bins)
        
    master_bins = pd.concat(all_bins).reset_index(drop=True)
    master_bins['bin_id'] = master_bins.index
    return master_bins

# --- 2. The Bridge Logic ---
def rasterize_labels(master_bins, bed_paths, cell_types, mode='hard'):
    """
    Overlaps BED intervals with Bins and resolves values.
    mode: 'hard' (Priority Max) or 'soft' (Weighted Average)
    """
    
    # Prepare Master Labels
    master_labels = master_bins[['bin_id', 'chrom', 'start', 'end']].copy()
    bin_gr = pr.PyRanges(master_labels.rename(columns={"start": "Start", "end": "End", "chrom": "Chromosome"}))
    
    # Priority for Hard Labels (Higher overrides lower)
    # 2=Exon, 1=Intron, 0=Desert
    PRIORITY = {2: 10, 1: 1, 0: 0} 

    for ct in cell_types:
        bed_path = bed_paths.get(ct)
        if not bed_path or not os.path.exists(bed_path):
            print(f"[WARN] Missing BED for {ct}, skipping...")
            continue
            
        print(f"[INFO] Processing {ct} ({mode})...")
        
        # Load BED
        # Headerless: chr, start, end, value
        df = pd.read_csv(bed_path, sep='\t', names=['Chromosome', 'Start', 'End', 'Value'])
        gr = pr.PyRanges(df)
        
        # Join Bins with Labels
        # This gives us the intersection of every bin and every label
        joined = bin_gr.join(gr).as_df()
        
        # Calculate overlap length for weighting
        joined['overlap_len'] = joined['End'] - joined['Start']
        
        if mode == 'hard':
            # Priority Resolution
            # If a bin hits multiple labels (edge case), take highest priority
            joined['prio'] = joined['Value'].map(PRIORITY).fillna(0)
            
            # Group by Bin -> Max Priority -> Get corresponding Value
            # (Sort by priority desc, take first)
            best_match = joined.sort_values('prio', ascending=False).drop_duplicates(subset=['bin_id'])
            final_series = best_match.set_index('bin_id')['Value']
            
        elif mode == 'soft':
            # Weighted Average Resolution
            # Sum(Value * Length) / Sum(Length)
            # Note: This assumes the label covers the overlap.
            
            joined['weighted_val'] = joined['Value'] * joined['overlap_len']
            
            grouped = joined.groupby('bin_id')
            total_val = grouped['weighted_val'].sum()
            total_len = grouped['overlap_len'].sum()
            
            avg_val = total_val / total_len
            final_series = avg_val
            
        # Merge into Master
        col_name = f"{ct}_Class"
        master_labels[col_name] = master_labels['bin_id'].map(final_series).fillna(0)
        
    return master_labels.drop(columns=['bin_id'])

# --- 3. Main Runner ---
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-gtf", required=True)
    parser.add_argument("--bin-size", type=int, default=32)
    parser.add_argument("--mode", choices=['hard', 'soft'], required=True)
    parser.add_argument("--out-tsv", required=True)
    
    # Dynamic arguments for cell types
    parser.add_argument("--b-bed", required=True)
    parser.add_argument("--t-bed", required=True)
    parser.add_argument("--m-bed", required=True)
    
    args = parser.parse_args()
    
    # Generate Bins
    master_bins = generate_bins(args.ref_gtf, args.bin_size)
    
    # Map Paths
    paths = {
        'B_Cell': args.b_bed,
        'T_Cell': args.t_bed,
        'Macrophage': args.m_bed
    }
    
    # Rasterize
    final_df = rasterize_labels(master_bins, paths, ['B_Cell', 'T_Cell', 'Macrophage'], mode=args.mode)
    
    # Save
    print(f"[INFO] Saving dense targets to {args.out_tsv}...")
    final_df.to_csv(args.out_tsv, sep='\t', index=False)

if __name__ == "__main__":
    main()