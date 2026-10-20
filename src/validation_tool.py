import argparse
import pandas as pd
import pyranges as pr
from generate_splicing_labels import load_stringtie_gtf

def parse_args():
    parser = argparse.ArgumentParser(description="Generate validation tracks comparing StringTie vs New Labels")
    parser.add_argument("--gtf", required=True, help="Original Dominant StringTie GTF")
    parser.add_argument("--hard-bed", required=True, help="The generated Hard Labels BED file")
    parser.add_argument("--out-bed", required=True, help="Output Validation BED file")
    return parser.parse_args()

def main():
    args = parse_args()
    
    print("[INFO] Loading Original StringTie GTF...")
    # Load original exons
    st_gr = load_stringtie_gtf(args.gtf)
    st_exons = st_gr[st_gr.Feature == 'exon']
    # Keep only coords
    st_exons = st_exons[['Chromosome', 'Start', 'End', 'Strand']]

    print("[INFO] Loading New Hard Labels...")
    # Load new labels (Label 2 = Exon)
    # BED file has no header: chrom, start, end, label
    new_df = pd.read_csv(args.hard_bed, sep='\t', names=['Chromosome', 'Start', 'End', 'label'])
    if 'Strand' not in new_df.columns:
        new_df['Strand'] = '.'
    
    # Filter for Exons (Label 2) in the new dataset
    new_exons = pr.PyRanges(new_df[new_df['label'] == 2])
    
    print("[INFO] Calculating Differences...")
    
    # 1. REMOVED EXONS (Red)
    # Present in StringTie, MISSING in New Labels
    removed_exons = st_exons.subtract(new_exons)
    if len(removed_exons) > 0:
        removed_exons.itemRgb = "255,0,0" # Red
        removed_exons.Name = "Removed_FP"
    
    # 2. ADDED EXONS (Green)
    # Present in New Labels, MISSING in StringTie
    added_exons = new_exons.subtract(st_exons)
    if len(added_exons) > 0:
        added_exons.itemRgb = "0,255,0" # Green
        added_exons.Name = "Added_FN"
        
    # Combine
    diff_gr = pr.concat([removed_exons, added_exons])
    
    if len(diff_gr) == 0:
        print("[WARN] No differences found! The datasets are identical.")
        return

    # Sort
    diff_gr = diff_gr.sort()
    
    # Output simple BED for IGV
    # BED9 Format for colors: chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb
    df_out = diff_gr.df.copy()
    
    # Fill BED fields
    df_out['score'] = 0
    if 'Strand' not in df_out.columns: df_out['Strand'] = '.'
    df_out['thickStart'] = df_out['Start']
    df_out['thickEnd'] = df_out['End']
    
    final_cols = ['Chromosome', 'Start', 'End', 'Name', 'score', 'Strand', 'thickStart', 'thickEnd', 'itemRgb']
    
    print(f"[INFO] Writing {len(df_out)} changes to {args.out_bed}")
    df_out[final_cols].to_csv(args.out_bed, sep='\t', header=False, index=False, quoting=3)

if __name__ == "__main__":
    main()