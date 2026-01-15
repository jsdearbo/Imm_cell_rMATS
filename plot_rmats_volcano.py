#!/usr/bin/env python3
"""Generate volcano plots for rMATS comparison outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from adjustText import adjust_text

    HAVE_ADJUST = True
except ImportError:  # pragma: no cover - optional dependency
    HAVE_ADJUST = False

EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")
DEFAULT_EXCLUDES = {"Arap1"}


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create volcano plots for each rMATS event type in comparison result folders."
    )
    parser.add_argument(
        "result_dirs",
        nargs="+",
        type=Path,
        help="One or more rMATS comparison directories (each containing *.MATS.JC.txt files).",
    )
    parser.add_argument(
        "--highlight",
        metavar="GENE",
        nargs="*",
        default=[],
        help="Genes to force-label on plots (repeat flag or space-separate).",
    )
    parser.add_argument(
        "--exclude",
        metavar="GENE",
        nargs="*",
        default=[],
        help="Genes to remove from plotting (repeat flag or space-separate).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of lowest-FDR events to annotate (default: %(default)s).",
    )
    parser.add_argument(
        "--output-subdir",
        default="plots",
        help="Folder created inside each comparison directory for plot outputs (default: %(default)s).",
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="Optional filename prefix prepended to each saved plot.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=160,
        help="Figure DPI for saved images (default: %(default)s).",
    )
    return parser.parse_args(argv)


def load_event_table(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None
    return df


def prep_volcano_df(df: pd.DataFrame, exclude_genes: set[str]) -> pd.DataFrame:
    working = df.assign(
        FDR=pd.to_numeric(df.get("FDR"), errors="coerce"),
        IncLevelDifference=pd.to_numeric(df.get("IncLevelDifference"), errors="coerce"),
    )
    filtered = working.dropna(subset=["FDR", "IncLevelDifference", "geneSymbol"])
    if exclude_genes:
        filtered = filtered[~filtered["geneSymbol"].isin(exclude_genes)]
    if filtered.empty:
        return filtered
    return filtered.assign(logFDR=lambda d: -np.log10(d["FDR"].clip(lower=1e-300)))


def annotate_points(ax, plot_df: pd.DataFrame, labels: pd.DataFrame) -> None:
    texts = []
    for _, row in labels.iterrows():
        ax.scatter(row["IncLevelDifference"], row["logFDR"], color="crimson", s=28, zorder=3)
        texts.append(
            ax.annotate(
                row["geneSymbol"],
                xy=(row["IncLevelDifference"], row["logFDR"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
                arrowprops=dict(arrowstyle="-", color="black", lw=0.4),
            )
        )
    if HAVE_ADJUST and texts:
        adjust_text(
            texts,
            ax=ax,
            force_points=0.05,
            force_text=0.05,
            expand_points=(1.02, 1.02),
            expand_text=(1.02, 1.02),
            only_move={"points": "y", "text": "xy"},
        )


def build_label_set(
    plot_df: pd.DataFrame,
    top_n: int,
    highlight_genes: set[str],
) -> pd.DataFrame:
    frames = []
    if top_n > 0 and not plot_df.empty:
        frames.append(plot_df.nsmallest(top_n, "FDR"))
    if highlight_genes:
        frames.append(plot_df[plot_df["geneSymbol"].isin(highlight_genes)])
    if not frames:
        return pd.DataFrame(columns=plot_df.columns)
    combined = pd.concat(frames, ignore_index=True)
    return combined.drop_duplicates(subset="geneSymbol")


def plot_event(
    comparison_dir: Path,
    event_type: str,
    df: pd.DataFrame,
    *,
    dpi: int,
    top_n: int,
    highlight_genes: set[str],
    exclude_genes: set[str],
    output_dir: Path,
    prefix: str,
) -> None:
    plot_df = prep_volcano_df(df, exclude_genes)
    if plot_df.empty:
        print(f"[skip] {comparison_dir.name} {event_type}: no valid rows after filtering", file=sys.stderr)
        return

    label_df = build_label_set(plot_df, top_n, highlight_genes)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(
        plot_df["IncLevelDifference"],
        plot_df["logFDR"],
        s=12,
        alpha=0.6,
        color="steelblue",
        edgecolor="none",
    )
    ax.axvline(0, color="gray", lw=1, ls="--")
    ax.set_xlabel("Δ PSI (IncLevelDifference)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(f"{event_type} events volcano ({comparison_dir.name.replace('_', ' ')})")

    if not label_df.empty:
        annotate_points(ax, plot_df, label_df)

    plt.tight_layout()
    filename = f"{prefix + '_' if prefix else ''}{comparison_dir.name}_{event_type}_volcano.png"
    output_path = output_dir / filename
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    print(f"[ok] wrote {output_path}")


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    exclude_genes = DEFAULT_EXCLUDES.union({gene.strip() for gene in args.exclude if gene.strip()})
    highlight_genes = {gene.strip() for gene in args.highlight if gene.strip()}

    for result_dir in args.result_dirs:
        if not result_dir.exists():
            print(f"[warn] {result_dir} missing; skipping", file=sys.stderr)
            continue
        output_dir = result_dir / args.output_subdir
        output_dir.mkdir(parents=True, exist_ok=True)
        for event in EVENT_TYPES:
            table_path = result_dir / f"{event}.MATS.JC.txt"
            table = load_event_table(table_path)
            if table is None:
                print(f"[skip] {table_path} not found or empty", file=sys.stderr)
                continue
            plot_event(
                result_dir,
                event,
                table,
                dpi=args.dpi,
                top_n=args.top_n,
                highlight_genes=highlight_genes,
                exclude_genes=exclude_genes,
                output_dir=output_dir,
                prefix=args.prefix,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
