import argparse
from pathlib import Path
from generate_splicing_labels import (
    load_rmats_se,
    load_rmats_ri,
    load_rmats_alt_site,
    load_stringtie_gtf,
    generate_hard_labels,
)

DEFAULTS = {
    "rmats_se": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_only/SE.MATS.JC.txt",
    "rmats_ri": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_only/RI.MATS.JC.txt",
    "rmats_a3ss": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_only/A3SS.MATS.JC.txt",
    "rmats_a5ss": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/rmats/B_only/A5SS.MATS.JC.txt",
    "gtf": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/compute_psi/stringtie_out/dominant_isoforms/B_Cell_dominant.gtf",
    "min_coverage": 20,
    "psi_low": 0.2,
    "psi_high": 0.8,
    "out_bed": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/compute_psi/hard_and_soft_label_testing/hard_output_test/B_cell_hard_labels.bed",
    "summary": "/bobross/jdearborn/public_data/ImmGen_Deep_RNAseq/compute_psi/hard_and_soft_label_testing/hard_output_test/B_cell_hard_label_summary.tsv",
}

def parse_args():
    parser = argparse.ArgumentParser(description="Test the hard-label correction flow.")
    parser.add_argument("--rmats-se", default=DEFAULTS["rmats_se"], help="Path to SE.MATS.JC.txt")
    parser.add_argument("--rmats-ri", default=DEFAULTS["rmats_ri"], help="Path to RI.MATS.JC.txt")
    parser.add_argument("--rmats-a3ss", default=DEFAULTS["rmats_a3ss"], help="Path to A3SS.MATS.JC.txt")
    parser.add_argument("--rmats-a5ss", default=DEFAULTS["rmats_a5ss"], help="Path to A5SS.MATS.JC.txt")
    parser.add_argument("--gtf", default=DEFAULTS["gtf"], help="Dominant StringTie GTF")
    parser.add_argument("--min-coverage", type=int, default=DEFAULTS["min_coverage"], help="rMATS coverage cutoff")
    parser.add_argument("--psi-low", type=float, default=DEFAULTS["psi_low"], help="Low PSI threshold")
    parser.add_argument("--psi-high", type=float, default=DEFAULTS["psi_high"], help="High PSI threshold")
    parser.add_argument("--out-bed", default=DEFAULTS["out_bed"], help="Output BED")
    parser.add_argument("--summary", default=DEFAULTS["summary"], help="Optional TSV summary")
    return parser.parse_args()

def main():
    args = parse_args()
    print(f"[INFO] Loading data...")
    rmats_se = load_rmats_se(args.rmats_se, args.min_coverage)
    rmats_ri = load_rmats_ri(args.rmats_ri, args.min_coverage)
    rmats_a3ss = load_rmats_alt_site(args.rmats_a3ss, args.min_coverage)
    rmats_a5ss = load_rmats_alt_site(args.rmats_a5ss, args.min_coverage)
    stringtie_gr = load_stringtie_gtf(args.gtf)

    print("[INFO] Generating Labels...")
    df = generate_hard_labels(
        stringtie_gr, 
        rmats_se, 
        rmats_ri, 
        rmats_a3ss, 
        rmats_a5ss,
        psi_low=args.psi_low, 
        psi_high=args.psi_high,
        min_fragment_len=5 # Added protection against subtraction artifacts
    )

    # Format as BED6: chr, start, end, name(label), score(0), strand
    # Ensuring label is string for Name column
    df['name'] = df['label'].astype(str)
    df['score'] = 0
    
    output_cols = ['Chromosome', 'Start', 'End', 'name', 'score', 'Strand']

    out_path = Path(args.out_bed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    df[output_cols].to_csv(out_path, sep="\t", header=False, index=False)

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        df['label'].value_counts().rename_axis("label").to_frame("count").to_csv(summary_path, sep="\t")

    print(f"[INFO] Wrote {len(df)} regions to {out_path}")

if __name__ == "__main__":
    main()