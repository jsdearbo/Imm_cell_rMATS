#!/usr/bin/env python3
"""
Generate volcano plots for rMATS comparison outputs.
Updated to prevent Y-axis squashing by extreme FDR outliers.
"""

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
except ImportError:
    HAVE_ADJUST = False

EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS", "MXE")
DEFAULT_EXCLUDES = {"Arap1"}

DEFAULT_Y_CAP = 15
COMBINED_EVENT_TYPES = ("SE", "RI", "A3SS", "A5SS")
DEFAULT_COMBINED_NAME = "all_events"
DEFAULT_TOP_CSV = 50

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
        default="plots_v2",
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
    parser.add_argument(
        "--y-cap",
        type=float,
        default=DEFAULT_Y_CAP,
        help="Cap -log10(FDR) for plotting (default: %(default)s).",
    )
    parser.add_argument(
        "--top-csv",
        type=int,
        default=DEFAULT_TOP_CSV,
        help="Top N rows to write to CSV per plot (default: %(default)s).",
    )
    parser.add_argument(
        "--combined",
        action="store_true",
        help="Also plot combined SE/RI/A3SS/A5SS events.",
    )
    parser.add_argument(
        "--combined-name",
        default=DEFAULT_COMBINED_NAME,
        help="Filename label for combined plot outputs (default: %(default)s).",
    )
    parser.add_argument(
        "--min-delta-psi",
        type=float,
        default=0.0,
        help="Minimum |IncLevelDifference| to include in plot (default: %(default)s).",
    )
    parser.add_argument(
        "--sig-fdr",
        type=float,
        default=0.05,
        help="FDR threshold for significance coloring (default: %(default)s).",
    )
    return parser.parse_args(argv)


def load_event_table(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None
    return df


def prep_volcano_df(
    df: pd.DataFrame, exclude_genes: set[str], min_delta_psi: float = 0.0
) -> pd.DataFrame:
    working = df.assign(
        FDR=pd.to_numeric(df.get("FDR"), errors="coerce"),
        IncLevelDifference=pd.to_numeric(df.get("IncLevelDifference"), errors="coerce"),
    )
    filtered = working.dropna(subset=["FDR", "IncLevelDifference", "geneSymbol"])

    if exclude_genes:
        filtered = filtered[~filtered["geneSymbol"].isin(exclude_genes)]

    if min_delta_psi > 0:
        filtered = filtered[filtered["IncLevelDifference"].abs() >= min_delta_psi]

    if filtered.empty:
        return filtered

    return filtered.assign(logFDR=lambda d: -np.log10(d["FDR"].clip(lower=1e-300)))


def annotate_points(
    ax,
    labels: pd.DataFrame,
    avoid_x: np.ndarray | None = None,
    avoid_y: np.ndarray | None = None,
) -> None:
    texts = []
    for _, row in labels.iterrows():
        x_val = row["IncLevelDifference"]
        y_val = row.get("logFDR_capped", row["logFDR"])
        is_capped = bool(row.get("is_capped", False))
        ax.scatter(
            x_val, y_val, color="crimson",
            s=35 if is_capped else 28,
            marker="v" if is_capped else "o",
            zorder=4,
        )
        texts.append(ax.text(x_val, y_val, row["geneSymbol"], fontsize=8))
    if HAVE_ADJUST and texts:
        adjust_kwargs: dict = dict(
            ax=ax,
            force_text=(0.5, 0.5),
            force_points=(0.5, 0.5),
            expand=(1.2, 1.4),
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5, shrinkA=5),
        )
        if avoid_x is not None and avoid_y is not None:
            adjust_kwargs["x"] = avoid_x
            adjust_kwargs["y"] = avoid_y
        adjust_text(texts, **adjust_kwargs)


def build_label_set(
    plot_df: pd.DataFrame,
    top_n: int,
    highlight_genes: set[str],
) -> pd.DataFrame:
    frames = []
    if top_n > 0 and not plot_df.empty:
        ranked = (
            plot_df.assign(_abs_delta=plot_df["IncLevelDifference"].abs())
            .sort_values(["FDR", "_abs_delta"], ascending=[True, False])
            .drop_duplicates(subset="geneSymbol")  # deduplicate first so top-N is always N unique genes
            .drop(columns="_abs_delta")
        )
        frames.append(ranked.head(top_n))
    if highlight_genes:
        frames.append(plot_df[plot_df["geneSymbol"].isin(highlight_genes)])

    if not frames:
        return pd.DataFrame(columns=plot_df.columns)

    combined = pd.concat(frames, ignore_index=True)
    return combined.drop_duplicates(subset="geneSymbol")


def build_combined_table(
    event_tables: dict[str, pd.DataFrame],
    event_types: tuple[str, ...],
) -> pd.DataFrame | None:
    frames = []
    for event in event_types:
        table = event_tables.get(event)
        if table is None:
            continue
        frames.append(table.assign(EventType=event))
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def write_top_csv(plot_df: pd.DataFrame, output_path: Path, top_n: int) -> None:
    if top_n <= 0 or plot_df.empty:
        return
    top_df = plot_df.nsmallest(top_n, "FDR").sort_values("FDR", ascending=True)
    top_df.to_csv(output_path, index=False)
    print(f"[ok] wrote {output_path}")


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
    y_cap: float | None,
    top_csv: int,
    min_delta_psi: float = 0.0,
    sig_fdr: float = 0.05,
    title: str | None = None,
    file_stem: str | None = None,
) -> None:
    plot_df = prep_volcano_df(df, exclude_genes, min_delta_psi)
    if plot_df.empty:
        print(f"[skip] {comparison_dir.name} {event_type}: no valid rows after filtering", file=sys.stderr)
        return

    label_df = build_label_set(plot_df, top_n, highlight_genes)

    if title is None:
        title = f"{event_type} events volcano ({comparison_dir.name.replace('_', ' ')})"
    if file_stem is None:
        base_prefix = f"{prefix}_" if prefix else ""
        file_stem = f"{base_prefix}{comparison_dir.name}_{event_type}_volcano"

    if y_cap is not None:
        plot_df = plot_df.assign(
            logFDR_capped=lambda d: np.minimum(d["logFDR"], y_cap),
            is_capped=lambda d: d["logFDR"] > y_cap,
        )
        label_df = label_df.assign(
            logFDR_capped=lambda d: np.minimum(d["logFDR"], y_cap),
            is_capped=lambda d: d["logFDR"] > y_cap,
        )
    else:
        plot_df = plot_df.assign(logFDR_capped=plot_df["logFDR"], is_capped=False)
        label_df = label_df.assign(logFDR_capped=label_df["logFDR"], is_capped=False)

    max_log_fdr = plot_df["logFDR"].max()

    fig, ax = plt.subplots(figsize=(6, 5))

    sig_mask = plot_df["FDR"] < sig_fdr
    non_sig = plot_df[~sig_mask]
    sig = plot_df[sig_mask]

    ax.scatter(
        non_sig["IncLevelDifference"], non_sig["logFDR_capped"],
        s=12, alpha=0.35, color="lightgray", edgecolor="none",
    )
    ax.scatter(
        sig["IncLevelDifference"], sig["logFDR_capped"],
        s=14, alpha=0.7, color="steelblue", edgecolor="none",
    )

    capped_df = plot_df[plot_df["is_capped"]]
    if not capped_df.empty:
        ax.scatter(
            capped_df["IncLevelDifference"],
            capped_df["logFDR_capped"],
            s=18,
            marker="v",
            color="black",
            alpha=0.7,
            label=f"Capped (n={len(capped_df)})",
        )
    
    ax.axvline(0, color="gray", lw=1, ls="--")

    if min_delta_psi > 0:
        for sign in (-1, 1):
            ax.axvline(sign * min_delta_psi, color="dimgray", lw=0.8, ls=":", zorder=1)
        ax.text(
            min_delta_psi, 0.03, f"|ΔPSI|≥{min_delta_psi}",
            transform=ax.get_xaxis_transform(),
            va="bottom", ha="left", fontsize=7, color="dimgray",
        )

    sig_y = -np.log10(sig_fdr)
    ax.axhline(sig_y, color="dimgray", lw=0.8, ls=":", zorder=1)
    ax.text(
        0.01, sig_y, f"FDR={sig_fdr}",
        transform=ax.get_yaxis_transform(),
        va="bottom", ha="left", fontsize=7, color="dimgray",
    )

    if y_cap is not None and max_log_fdr > y_cap:
        n_capped_labels = int(label_df["is_capped"].sum()) if not label_df.empty else 0
        label_headroom = max(2.0, 0.4 * max(n_capped_labels, 4))
        ax.set_ylim(bottom=-0.5, top=y_cap + label_headroom)
        # Remove ticks in the headroom region; mark the cap tick as ≥{y_cap}
        step = 5 if y_cap >= 10 else 2
        ticks = list(np.arange(0, y_cap + 0.001, step))
        ax.set_yticks(ticks)
        ax.set_yticklabels([
            f"≥{int(t)}" if t == y_cap else str(int(t)) for t in ticks
        ])

    ax.set_xlabel("Δ PSI (IncLevelDifference)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(title, pad=12)

    parts = comparison_dir.name.split("_vs_")
    if len(parts) == 2:
        group1, group2 = (p.replace("_", " ") for p in parts)
        ax.text(
            0.02, 0.03, f"← Higher in {group2}",
            transform=ax.transAxes, fontsize=7.5, ha="left", va="bottom",
            color="dimgray", style="italic",
        )
        ax.text(
            0.98, 0.03, f"Higher in {group1} →",
            transform=ax.transAxes, fontsize=7.5, ha="right", va="bottom",
            color="dimgray", style="italic",
        )

    if not label_df.empty:
        annotate_points(
            ax, label_df,
            avoid_x=plot_df["IncLevelDifference"].values,
            avoid_y=plot_df["logFDR_capped"].values,
        )

    if not capped_df.empty:
        ax.legend(loc="upper right", frameon=False, fontsize=8)

    plt.tight_layout()
    output_path = output_dir / f"{file_stem}.png"
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    print(f"[ok] wrote {output_path}")

    csv_path = output_dir / f"{file_stem}_top{top_csv}.csv"
    write_top_csv(plot_df, csv_path, top_csv)


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
        event_tables: dict[str, pd.DataFrame] = {}
        for event in EVENT_TYPES:
            table_path = result_dir / f"{event}.MATS.JC.txt"
            table = load_event_table(table_path)
            if table is None:
                print(f"[skip] {table_path} not found or empty", file=sys.stderr)
                continue
            event_tables[event] = table
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
                y_cap=args.y_cap,
                top_csv=args.top_csv,
                min_delta_psi=args.min_delta_psi,
                sig_fdr=args.sig_fdr,
            )
        combined_table = build_combined_table(event_tables, COMBINED_EVENT_TYPES)
        if combined_table is None:
            print(
                f"[skip] {result_dir} combined: no valid rows after filtering",
                file=sys.stderr,
            )
        else:
            base_prefix = f"{args.prefix}_" if args.prefix else ""
            combined_stem = f"{base_prefix}{result_dir.name}_{args.combined_name}_volcano"
            combined_title = (
                "All events (SE, RI, A3SS, A5SS) volcano "
                f"({result_dir.name.replace('_', ' ')})"
            )
            plot_event(
                result_dir,
                args.combined_name,
                combined_table,
                dpi=args.dpi,
                top_n=args.top_n,
                highlight_genes=highlight_genes,
                exclude_genes=exclude_genes,
                output_dir=output_dir,
                prefix=args.prefix,
                y_cap=args.y_cap,
                top_csv=args.top_csv,
                min_delta_psi=args.min_delta_psi,
                sig_fdr=args.sig_fdr,
                title=combined_title,
                file_stem=combined_stem,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))