#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import re


def make_manifest(
    soft_label_root: str,
    out_tsv: str,
    bed_suffix: str = ".bed",
):
    """
    Build manifest.tsv listing all per-sample soft-label BEDs.

    Output columns:
      sample_id, cell_type, bed_path
    """

    soft_label_root = Path(soft_label_root)
    rows = []

    cell_dirs = sorted(p for p in soft_label_root.iterdir() if p.is_dir())

    if cell_dirs:
        for cell_dir in cell_dirs:
            cell_type = cell_dir.name

            bed_files = sorted(cell_dir.glob(f"*{bed_suffix}"))
            if not bed_files:
                print(f"[WARN] No BED files found in {cell_dir}")
                continue

            for bed_path in bed_files:
                name = bed_path.stem

                # Try to extract replicate/sample index if present
                m = re.search(r"psi_(\d+)", name)
                replicate_id = m.group(1) if m else None

                sample_id = f"{cell_type}_{name}"

                rows.append({
                    "sample_id": sample_id,
                    "cell_type": cell_type,
                    "bed_path": str(bed_path.resolve()),
                    "replicate_id": replicate_id,
                })
    else:
        bed_files = sorted(soft_label_root.glob(f"*{bed_suffix}"))
        if not bed_files:
            raise RuntimeError(f"No BED files found in {soft_label_root}")

        for bed_path in bed_files:
            name = bed_path.stem

            m = re.search(r"psi_(\d+)", name)
            replicate_id = m.group(1) if m else None

            cell_type = name.split("_soft", 1)[0]
            sample_id = f"{cell_type}_{name}"

            rows.append({
                "sample_id": sample_id,
                "cell_type": cell_type,
                "bed_path": str(bed_path.resolve()),
                "replicate_id": replicate_id,
            })

    if not rows:
        raise RuntimeError("No samples found — manifest would be empty")

    df = pd.DataFrame(rows)

    # Stable ordering = reproducibility
    df = df.sort_values(
        ["cell_type", "replicate_id", "sample_id"],
        na_position="last"
    ).reset_index(drop=True)

    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(out_tsv, sep="\t", index=False)

    print(f"[INFO] Wrote manifest with {len(df)} samples → {out_tsv}")


def parse_args():
    ap = argparse.ArgumentParser(description="Create manifest.tsv from per-sample soft-label BEDs.")
    ap.add_argument("--soft-label-root", required=True,
                    help="Root directory containing per-cell-type BED subdirectories")
    ap.add_argument("--out", required=True,
                    help="Output manifest.tsv path")
    return ap.parse_args()


if __name__ == "__main__":
    args = parse_args()
    make_manifest(
        soft_label_root=args.soft_label_root,
        out_tsv=args.out,
    )
