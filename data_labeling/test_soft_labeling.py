import argparse
from pathlib import Path
from generate_splicing_labels import (
    load_rmats_se,
    load_rmats_ri,
    load_rmats_alt_site,
    load_stringtie_gtf,
    generate_soft_labels,
)

DEFAULTS = {
    "rmats_se": "/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/results/B_only/SE.MATS.JC.txt",
    "rmats_ri": "/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/results/B_only/RI.MATS.JC.txt",
    "rmats_a3ss": "/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/results/B_only/A3SS.MATS.JC.txt",
    "rmats_a5ss": "/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/results/B_only/A5SS.MATS.JC.txt",
    "gtf": "/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/compute_psi/stringtie_out/dominant_isoforms/B_Cell_dominant.gtf",
    "min_coverage": 20,
    "out_bed": "/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/test_soft_labeling_output/B_cell_soft_labels.bed",
    "summary": "/mnt/data_sda/jdearborn/scratch/Imm_cell_rMATS/data_labeling/test_soft_labeling_output/B_cell_soft_label_summary.tsv",
}

def parse_args():
    parser = argparse.ArgumentParser(description="Test the soft-label correction flow.")
    parser.add_argument("--rmats-se", default=DEFAULTS["rmats_se"], help="Path to SE.MATS.JC.txt")
    parser.add_argument("--rmats-ri", default=DEFAULTS["rmats_ri"], help="Path to RI.MATS.JC.txt")
    parser.add_argument("--rmats-a3ss", default=DEFAULTS["rmats_a3ss"], help="Path to A3SS.MATS.JC.txt")
    parser.add_argument("--rmats-a5ss", default=DEFAULTS["rmats_a5ss"], help="Path to A5SS.MATS.JC.txt")
    parser.add_argument("--gtf", default=DEFAULTS["gtf"], help="Dominant StringTie GTF")
    parser.add_argument("--min-coverage", type=int, default=DEFAULTS["min_coverage"], help="rMATS coverage cutoff")
    parser.add_argument("--out-bed", default=DEFAULTS["out_bed"], help="Output BED (chr start end label)")
    parser.add_argument("--summary", default=DEFAULTS["summary"], help="Optional TSV summary of counts by label")
    return parser.parse_args()

def main():
    args = parse_args()
    print(f"[INFO] Loading SE from {args.rmats_se}")
    rmats_se = load_rmats_se(args.rmats_se, args.min_coverage)
    
    print(f"[INFO] Loading RI from {args.rmats_ri}")
    rmats_ri = load_rmats_ri(args.rmats_ri, args.min_coverage)

    print(f"[INFO] Loading A3SS from {args.rmats_a3ss}")
    rmats_a3ss = load_rmats_alt_site(args.rmats_a3ss, args.min_coverage)

    print(f"[INFO] Loading A5SS from {args.rmats_a5ss}")
    rmats_a5ss = load_rmats_alt_site(args.rmats_a5ss, args.min_coverage)
    
    print(f"[INFO] Loading GTF from {args.gtf}")
    stringtie_gr = load_stringtie_gtf(args.gtf)

    print("[INFO] Generating Labels...")
    df = generate_soft_labels(
        stringtie_gr,
        rmats_se,
        rmats_ri=rmats_ri,
        rmats_a3ss=rmats_a3ss,
        rmats_a5ss=rmats_a5ss,
    )

    out_path = Path(args.out_bed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Add default confidence if not present
    if 'confidence' not in df.columns:
        df['confidence'] = 0.0
    
    df[['Chromosome', 'Start', 'End', 'p_exon', 'confidence']].to_csv(
        out_path, sep="\t", header=False, index=False
    )

    if args.summary:
        summary_path = Path(args.summary)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        df['p_exon'].describe().rename_axis('stat').to_frame('value').to_csv(
            summary_path, sep="\t"
        )

    print(f"[INFO] Wrote {len(df)} regions to {out_path}")

if __name__ == "__main__":
    main()