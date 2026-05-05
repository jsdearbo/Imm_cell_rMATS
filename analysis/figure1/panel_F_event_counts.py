"""
panel_F_event_counts.py

Figure 1F: rMATS alternative splicing event counts by cell type and event class.
Grouped bar chart: x = event type, hue = cell type.
Two bar series per group: total events and FDR < 0.05 events.

Output: figure1_out/panel_F_event_counts.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ─────────────────────────────────────────────────────────────────────
RMATS_BASE = os.path.join(os.path.dirname(__file__), "../../results")
OUT_DIR    = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

CELL_TYPES = {
    "B_Cell":     "B_only",
    "T_Cell":     "T_only",
    "Macrophage": "Macrophage_only",
}
EVENT_TYPES = ["SE", "RI", "A3SS", "A5SS", "MXE"]
MIN_READS   = 20   # minimum total junction reads (IJC + SJC) per event

COLORS = {
    "B_Cell":     "#2166AC",
    "T_Cell":     "#D6604D",
    "Macrophage": "#4DAC26",
}


def _parse_read_sum(series) -> pd.Series:
    """Sum comma-separated read counts per row."""
    def _sum(s):
        if pd.isna(s):
            return 0.0
        try:
            return sum(float(x) for x in str(s).split(",") if x not in ("", "NA"))
        except Exception:
            return 0.0
    return series.apply(_sum)


def count_events(cell_type, subdir, min_reads=MIN_READS):
    rows = []
    for event in EVENT_TYPES:
        path = os.path.join(RMATS_BASE, subdir, f"{event}.MATS.JC.txt")
        if not os.path.exists(path):
            print(f"  [WARN] missing: {path}")
            rows.append({"cell_type": cell_type, "event_type": event,
                         "total": 0, "covered": 0})
            continue
        df = pd.read_csv(path, sep="\t")
        total = len(df)

        # Coverage filter: sum of IJC_SAMPLE_1 + SJC_SAMPLE_1 across replicates
        ijc = _parse_read_sum(df["IJC_SAMPLE_1"]) if "IJC_SAMPLE_1" in df.columns else pd.Series(0.0, index=df.index)
        sjc = _parse_read_sum(df["SJC_SAMPLE_1"]) if "SJC_SAMPLE_1" in df.columns else pd.Series(0.0, index=df.index)
        covered = int(((ijc + sjc) >= min_reads).sum())

        rows.append({"cell_type": cell_type, "event_type": event,
                     "total": total, "covered": covered})
        print(f"  [{cell_type}] {event}: {total} total, {covered} with ≥{min_reads} reads")
    return rows


def plot(counts_df):
    ct_order    = ["B_Cell", "T_Cell", "Macrophage"]
    event_order = EVENT_TYPES

    n_events = len(event_order)
    n_ct     = len(ct_order)
    bar_w    = 0.22
    group_w  = n_ct * bar_w + 0.08

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.5), sharey=False)

    for ax_idx, (col, ylabel, title_suffix) in enumerate([
        ("total",   "Event count (×10³)", "all detected"),
        ("covered", "Event count (×10³)", f"≥{MIN_READS} junction reads"),
    ]):
        ax = axes[ax_idx]
        for ci, ct in enumerate(ct_order):
            sub = counts_df[counts_df["cell_type"] == ct]
            x_base = np.arange(n_events) * group_w
            x = x_base + ci * bar_w
            vals = [sub[sub["event_type"] == ev][col].values[0] / 1e3
                    if len(sub[sub["event_type"] == ev]) else 0
                    for ev in event_order]
            ax.bar(x, vals, width=bar_w, color=COLORS[ct], alpha=0.85,
                   label=ct.replace("_", " "), zorder=2)

        ax.set_xticks(np.arange(n_events) * group_w + bar_w)
        ax.set_xticklabels(event_order, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_title(f"rMATS events — {title_suffix}", fontsize=9)
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.grid(True, lw=0.5, alpha=0.4, zorder=0)
        ax.set_axisbelow(True)
        if ax_idx == 0:
            ax.legend(fontsize=7, framealpha=0.9, edgecolor="#CCCCCC")

    plt.tight_layout()
    out = os.path.join(OUT_DIR, "panel_F_event_counts.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    all_rows = []
    for ct, subdir in CELL_TYPES.items():
        print(f"Loading {ct}...")
        all_rows.extend(count_events(ct, subdir))

    counts_df = pd.DataFrame(all_rows)
    print("\n" + counts_df.pivot_table(
        index="cell_type", columns="event_type",
        values=["total", "covered"], aggfunc="sum"
    ).to_string())

    plot(counts_df)
