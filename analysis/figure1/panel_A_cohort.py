"""
panel_A_cohort.py

Figure 1A: RNA-seq cohort overview.
Horizontal Cleveland dot plot showing primary mapped read pairs (millions) per sample,
grouped by cell type, with excluded (non-PE150) samples shown in gray.

Output: figure1_out/panel_A_cohort.pdf
"""

import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ─────────────────────────────────────────────────────────────────────
FLAGSTAT_ROOT = os.path.join(
    os.path.dirname(__file__), "../../results/flagstat"
)
OUT_DIR = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Inclusion lists (from config/{B,T,Mac}_Cell_bams.txt) ────────────────────
INCLUDED = {
    "B_Cell":     {"SRR7988747", "SRR7988759", "SRR7988771", "SRR14714583", "SRR14714587"},
    "T_Cell":     {"SRR17143398", "SRR17143399", "SRR25152365", "SRR25152366"},
    "Macrophage": {"SRR14049271", "SRR14049272", "SRR14049273",
                   "SRR15402174", "SRR15402175", "SRR31668386", "SRR31668387"},
}

CELL_TYPE_DIRS = {
    "B_Cell":     "flagstat_results_b_cell",
    "T_Cell":     "flagstat_results_t_cell",
    "Macrophage": "flagstat_results_macrophage",
}

# Publication colors
COLORS = {
    "B_Cell":     "#2166AC",
    "T_Cell":     "#D6604D",
    "Macrophage": "#4DAC26",
    "excluded":   "#BDBDBD",
}


def parse_flagstat(path):
    """Return (primary_mapped, properly_paired_pct) from a flagstat file."""
    primary_mapped = None
    paired_pct = None
    with open(path) as fh:
        for line in fh:
            m = re.match(r"(\d+) \+ \d+ primary mapped", line)
            if m:
                primary_mapped = int(m.group(1))
            m2 = re.match(r"(\d+) \+ \d+ properly paired \((\d+\.\d+)%", line)
            if m2:
                paired_pct = float(m2.group(2))
    return primary_mapped, paired_pct


def load_all_samples():
    records = []
    for cell_type, subdir in CELL_TYPE_DIRS.items():
        d = os.path.join(FLAGSTAT_ROOT, subdir)
        for fname in sorted(os.listdir(d)):
            if not fname.endswith("_flagstat.txt"):
                continue
            srr = fname.replace("_flagstat.txt", "")
            primary_mapped, paired_pct = parse_flagstat(os.path.join(d, fname))
            included = srr in INCLUDED[cell_type]
            records.append({
                "srr":           srr,
                "cell_type":     cell_type,
                "primary_mapped": primary_mapped,
                "paired_pct":    paired_pct,
                "included":      included,
            })
    return records


def plot(records):
    # Sort: cell type order, then included-first within each type, then by read depth desc
    ct_order = ["B_Cell", "T_Cell", "Macrophage"]
    records.sort(key=lambda r: (
        ct_order.index(r["cell_type"]),
        0 if r["included"] else 1,
        -r["primary_mapped"],
    ))

    n = len(records)
    mapped_M = [r["primary_mapped"] / 1e6 for r in records]
    labels = [r["srr"] for r in records]
    colors = [COLORS[r["cell_type"]] if r["included"] else COLORS["excluded"]
              for r in records]
    alphas = [1.0 if r["included"] else 0.45 for r in records]

    fig, ax = plt.subplots(figsize=(5.5, 5.2))

    y_positions = np.arange(n)[::-1]  # top = first record

    for i, (y, m, c, a, rec) in enumerate(zip(y_positions, mapped_M, colors, alphas, records)):
        ax.plot([0, m], [y, y], color=c, alpha=a, lw=1.2, zorder=1)
        ax.scatter([m], [y], color=c, alpha=a, s=40, zorder=2,
                   marker="o" if rec["included"] else "x")

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Primary mapped read pairs (millions)", fontsize=8)
    ax.set_xlim(left=0)
    ax.tick_params(axis="x", labelsize=7)

    # Cell-type separators
    y_sep_positions = []
    prev_ct = None
    for i, rec in enumerate(records):
        if prev_ct is not None and rec["cell_type"] != prev_ct:
            y_sep_positions.append(n - i - 0.5)
        prev_ct = rec["cell_type"]
    for yp in y_sep_positions:
        ax.axhline(yp, color="#CCCCCC", lw=0.8, linestyle="--", zorder=0)

    # Cell-type labels on right margin
    ct_midpoints = {}
    for ct in ct_order:
        ys = [y_positions[i] for i, r in enumerate(records) if r["cell_type"] == ct]
        if ys:
            ct_midpoints[ct] = (min(ys) + max(ys)) / 2
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(list(ct_midpoints.values()))
    ax2.set_yticklabels(
        [ct.replace("_", " ") for ct in ct_order if ct in ct_midpoints],
        fontsize=8, fontweight="bold",
    )
    for tick, ct in zip(ax2.get_yticklabels(),
                        [c for c in ct_order if c in ct_midpoints]):
        tick.set_color(COLORS[ct])
    ax2.tick_params(length=0)

    # Legend
    legend_handles = [
        mpatches.Patch(color=COLORS["B_Cell"],     label="B Cell (n=5)"),
        mpatches.Patch(color=COLORS["T_Cell"],     label="T Cell (n=4)"),
        mpatches.Patch(color=COLORS["Macrophage"], label="Macrophage (n=7)"),
        mpatches.Patch(color=COLORS["excluded"],   label="Excluded (not PE150)",
                       alpha=0.5),
    ]
    ax.legend(handles=legend_handles, fontsize=6.5, loc="lower right",
              framealpha=0.9, edgecolor="#CCCCCC")

    ax.set_title("RNA-seq cohort (PE150, mm10)", fontsize=9, pad=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax2.spines[["top", "left"]].set_visible(False)

    plt.tight_layout()
    out = os.path.join(OUT_DIR, "panel_A_cohort.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    print(f"Saved: {out}")

    # Summary
    included = [r for r in records if r["included"]]
    print(f"\nIncluded samples: {len(included)}")
    for ct in ct_order:
        ct_inc = [r for r in included if r["cell_type"] == ct]
        reads_M = [r["primary_mapped"] / 1e6 for r in ct_inc]
        print(f"  {ct}: n={len(ct_inc)}, "
              f"mapped reads {min(reads_M):.0f}–{max(reads_M):.0f} M "
              f"(median {np.median(reads_M):.0f} M)")


if __name__ == "__main__":
    records = load_all_samples()
    plot(records)
