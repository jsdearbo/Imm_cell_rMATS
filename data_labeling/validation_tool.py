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
    st_gr = load_stringtie_gtf(args.gtf)
    st_exons = st_gr[st_gr.Feature == 'exon']
    st_exons = st_exons[['Chromosome', 'Start', 'End', 'Strand']]

    print("[INFO] Loading New Hard Labels...")
    # Updated to read standard BED6 format
    # Columns: Chromosome, Start, End, Name(label), Score, Strand
    try:
        new_df = pd.read_csv(
            args.hard_bed, 
            sep='\t', 
            names=['Chromosome', 'Start', 'End', 'label', 'score', 'Strand']
        )
    except Exception as e:
        print(f"[ERROR] Could not read BED file (expecting BED6 format): {e}")
        return

    # Filter for Exons (Label 2)
    # Note: label comes in as 'Name' column (string), so check for "2" or "2.0" or integer 2
    new_df['label'] = new_df['label'].astype(str)
    new_exons = pr.PyRanges(new_df[new_df['label'].str.startswith('2')])
    
    print(f"[INFO] Comparing {len(st_exons)} old exons vs {len(new_exons)} new exons...")
    
    # 1. REMOVED EXONS (Red)
    # Note: PyRanges subtract is strand-specific by default.
    removed_exons = st_exons.subtract(new_exons)
    if len(removed_exons) > 0:
        removed_exons.itemRgb = "255,0,0" # Red
        removed_exons.Name = "Removed_FP"
    
    # 2. ADDED EXONS (Green)
    added_exons = new_exons.subtract(st_exons)
    if len(added_exons) > 0:
        added_exons.itemRgb = "0,255,0" # Green
        added_exons.Name = "Added_FN"
        
    # Combine
    diff_gr = pr.concat([removed_exons, added_exons])
    
    if len(diff_gr) == 0:
        print("[WARN] No differences found! The datasets are identical.")
        return

    # Output BED9
    df_out = diff_gr.df.copy()
    df_out['score'] = 0
    df_out['thickStart'] = df_out['Start']
    df_out['thickEnd'] = df_out['End']
    
    final_cols = ['Chromosome', 'Start', 'End', 'Name', 'score', 'Strand', 'thickStart', 'thickEnd', 'itemRgb']
    
    print(f"[INFO] Writing {len(df_out)} changes to {args.out_bed}")
    df_out[final_cols].to_csv(args.out_bed, sep='\t', header=False, index=False, quoting=3)

if __name__ == "__main__":
    main()