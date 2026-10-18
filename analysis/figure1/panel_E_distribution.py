"""
panel_E_distribution.py

Figure 1E: Genome-wide PSI label distribution across cell types.
Two-part figure:
  Left (3 panels): Per-cell-type PSI histogram on log y-scale — reveals trimodal
                   structure (skipped / regulated / included).
  Right (1 panel): Stacked bar chart of category fractions per cell type.

Only bins where mask == 1 are counted (bin-level validity flag, not has_signal).

Output: figure1_out/panel_E_distribution.pdf
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import figstyle

# ── Paths ─────────────────────────────────────────────────────────────────────
SHARD_DIR = os.path.join(
    os.path.dirname(__file__),
    "../../data_labeling/results/labelled_data/multitask_balanced/train",
)
OUT_DIR = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

# Task ordering from NPZ celltype_idx_N keys: 0=B_Cell, 1=Macrophage, 2=T_Cell
CELL_TYPES = ["B Cell", "Macrophage", "T Cell"]
TASK_COLORS = {
    "B Cell":     "#2166AC",
    "Macrophage": "#4DAC26",
    "T Cell":     "#D6604D",
}

CAT_COLORS = {
    "skipped":   "#5E9ACC",   # PSI < 0.1
    "regulated": "#E67E22",   # PSI 0.1–0.9
    "included":  "#C0392B",   # PSI > 0.9
}

LO = 0.1
HI = 0.9
MAX_SHARDS = 10


def load_psi_values(shard_dir, max_shards=MAX_SHARDS):
    """Return dict cell_type → 1-D array of all masked PSI values."""
    shards = sorted(glob.glob(os.path.join(shard_dir, "*.npz")))[:max_shards]
    psi_per_ct = {ct: [] for ct in CELL_TYPES}
    for path in shards:
        d = np.load(path, allow_pickle=True)
        targets = d["targets"]   # (N, bins, 3)
        mask    = d["mask"]      # (N, bins, 3)
        for ti, ct in enumerate(CELL_TYPES):
            vals = targets[:, :, ti][mask[:, :, ti] == 1]
            psi_per_ct[ct].append(vals)
        print(f"  loaded {os.path.basename(path)}")
    return {ct: np.concatenate(arrs) for ct, arrs in psi_per_ct.items()}


def category_fractions(vals):
    lo  = (vals < LO).mean() * 100
    mid = ((vals >= LO) & (vals <= HI)).mean() * 100
    hi  = (vals > HI).mean() * 100
    return lo, mid, hi


def plot(psi_per_ct):
    figstyle.apply(compact=True)
    fig, axes = plt.subplots(
        1, 4,
        figsize=(10.5, 3.4),
        gridspec_kw={"width_ratios": [1, 1, 1, 0.85]},
    )
    fig.subplots_adjust(wspace=0.38)

    # PSI values cluster on multiples of 1/32 (label-averaging resolution), so
    # 60 arbitrary equal-width bins alias against that grid: some bins land on
    # a mass point, neighboring bins fall in the gap, producing an alternating
    # high/low "zebra stripe" artifact. Centering 32 bins on k/32 removes it.
    bins        = np.clip((np.arange(33) - 0.5) / 32, 0, 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_w       = np.diff(bins)

    cat_data = {}   # collected for the stacked bar panel

    # ── Histogram panels (log scale) ─────────────────────────────────────────
    for ax, ct in zip(axes[:3], CELL_TYPES):
        vals  = psi_per_ct[ct]
        color = TASK_COLORS[ct]

        counts, _ = np.histogram(vals, bins=bins, density=False)
        freq = counts / counts.sum()   # fraction per bin

        # Color bars by PSI category
        bar_colors = np.where(
            bin_centers < LO,  CAT_COLORS["skipped"],
            np.where(bin_centers > HI, CAT_COLORS["included"],
                     CAT_COLORS["regulated"])
        )
        ax.bar(bin_centers, freq, width=bin_w,
               color=bar_colors, alpha=0.85, lw=0, zorder=2)

        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_ylim(1e-5, 1)
        ax.set_xlabel("PSI")
        ax.set_title(ct, color=color, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.grid(True, lw=0.4, alpha=0.35, zorder=0, which="both")
        ax.set_axisbelow(True)

        lo_f, mid_f, hi_f = category_fractions(vals)
        cat_data[ct] = (lo_f, mid_f, hi_f)

    axes[0].set_ylabel("Fraction of labeled bins")
    axes[1].set_yticklabels([])
    axes[2].set_yticklabels([])

    # ── Stacked bar panel ─────────────────────────────────────────────────────
    ax = axes[3]
    ct_labels = CELL_TYPES
    x = np.arange(len(ct_labels))

    bottoms = np.zeros(len(ct_labels))
    for cat_key, cat_label in [
        ("skipped",   f"PSI < {LO}"),
        ("regulated", f"PSI {LO}–{HI}"),
        ("included",  f"PSI > {HI}"),
    ]:
        vals_cat = np.array([cat_data[ct][["skipped","regulated","included"].index(cat_key)]
                             for ct in ct_labels])
        bars = ax.bar(x, vals_cat, bottom=bottoms,
                      color=CAT_COLORS[cat_key], width=0.55,
                      label=cat_label, zorder=2)
        # Label segments that are large enough to read
        for xi, (bot, v) in enumerate(zip(bottoms, vals_cat)):
            if v >= 3.0:
                ax.text(xi, bot + v / 2, f"{v:.1f}%",
                        ha="center", va="center", fontsize=figstyle.ANNOT_XS,
                        color="white", fontweight="bold")
        bottoms += vals_cat

    ax.set_xticks(x)
    ax.set_xticklabels([ct.replace(" ", "\n") for ct in ct_labels])
    ax.set_ylabel("% of labeled bins")
    ax.set_ylim(0, 102)
    ax.set_title("PSI breakdown", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, lw=0.4, alpha=0.35, zorder=0)
    ax.set_axisbelow(True)

    # Shared legend (placed below stacked bar)
    legend_handles = [
        Patch(color=CAT_COLORS["skipped"],   label=f"PSI < {LO}  (skipped)"),
        Patch(color=CAT_COLORS["regulated"], label=f"PSI {LO}–{HI}  (regulated)"),
        Patch(color=CAT_COLORS["included"],  label=f"PSI > {HI}  (included)"),
    ]
    fig.legend(handles=legend_handles,
               loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.08),
               framealpha=0.9, edgecolor="#CCCCCC")

    fig.suptitle("Genome-wide PSI label distribution", y=1.02)

    out = os.path.join(OUT_DIR, "panel_E_distribution.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    print(f"Saved: {out}")

    print("\n=== PSI category breakdown ===")
    for ct in CELL_TYPES:
        vals = psi_per_ct[ct]
        lo, mid, hi = category_fractions(vals)
        print(f"  {ct}: {len(vals):,} masked bins | "
              f"<{LO}: {lo:.1f}% | {LO}–{HI}: {mid:.1f}% | >{HI}: {hi:.1f}%")


if __name__ == "__main__":
    print(f"Loading up to {MAX_SHARDS} train shards from:\n  {SHARD_DIR}\n")
    psi_per_ct = load_psi_values(SHARD_DIR)
    plot(psi_per_ct)
