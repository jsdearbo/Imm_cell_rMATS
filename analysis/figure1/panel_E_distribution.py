"""
panel_E_distribution.py

Figure 1E: Genome-wide PSI label distribution across cell types.
Histogram + KDE showing trimodal structure (skipped / regulated / included).
Only bins where mask == 1 are counted.

Loads a sample of train shards to keep runtime manageable.

Output: figure1_out/panel_E_distribution.pdf
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# ── Paths ─────────────────────────────────────────────────────────────────────
SHARD_DIR = os.path.join(
    os.path.dirname(__file__),
    "../../data_labeling/results/labelled_data/multitask_balanced/train",
)
OUT_DIR = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

# Task ordering from NPZ celltype_idx_N keys
CELL_TYPES = ["B Cell", "Macrophage", "T Cell"]  # tasks 0, 1, 2
TASK_COLORS = {
    "B Cell":     "#2166AC",
    "Macrophage": "#4DAC26",
    "T Cell":     "#D6604D",
}

# PSI category thresholds
LO = 0.1
HI = 0.9
MAX_SHARDS = 10   # sample shards for speed; increase for final figure


def load_psi_values(shard_dir, max_shards=MAX_SHARDS):
    """Return dict cell_type → 1-D array of all masked PSI values."""
    shards = sorted(glob.glob(os.path.join(shard_dir, "*.npz")))[:max_shards]
    psi_per_ct = {ct: [] for ct in CELL_TYPES}

    for path in shards:
        d = np.load(path, allow_pickle=True)
        targets = d["targets"]   # (N, bins, 3)
        mask    = d["mask"]      # (N, bins, 3)
        for ti, ct in enumerate(CELL_TYPES):
            vals = targets[:, :, ti][mask[:, :, ti].astype(bool)]
            psi_per_ct[ct].append(vals)
        print(f"  loaded {os.path.basename(path)}")

    return {ct: np.concatenate(arrs) for ct, arrs in psi_per_ct.items()}


def category_fractions(vals):
    lo = (vals < LO).mean() * 100
    mid = ((vals >= LO) & (vals <= HI)).mean() * 100
    hi  = (vals > HI).mean() * 100
    return lo, mid, hi


def plot(psi_per_ct):
    n_ct = len(CELL_TYPES)
    fig, axes = plt.subplots(1, n_ct, figsize=(8, 2.8), sharey=False)

    bins = np.linspace(0, 1, 60)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    for ax, ct in zip(axes, CELL_TYPES):
        vals = psi_per_ct[ct]
        color = TASK_COLORS[ct]

        counts, _ = np.histogram(vals, bins=bins, density=False)
        freq = counts / counts.sum()

        # Color bins by PSI category
        cat_colors = np.where(
            bin_centers < LO, "#5E9ACC",
            np.where(bin_centers > HI, "#C0392B", "#E67E22")
        )

        ax.bar(bin_centers, freq, width=bins[1] - bins[0],
               color=cat_colors, alpha=0.8, zorder=2)

        # KDE overlay
        if len(vals) > 100:
            kde = gaussian_kde(vals, bw_method=0.05)
            x_kde = np.linspace(0, 1, 400)
            y_kde = kde(x_kde)
            y_kde_scaled = y_kde * (bins[1] - bins[0])
            ax.plot(x_kde, y_kde_scaled, color=color, lw=1.5, zorder=3)

        lo_frac, mid_frac, hi_frac = category_fractions(vals)
        ax.set_title(ct, fontsize=8.5, color=color, fontweight="bold")
        ax.text(0.05, 0.97,
                f"<0.1: {lo_frac:.1f}%\n0.1–0.9: {mid_frac:.1f}%\n>0.9: {hi_frac:.1f}%",
                transform=ax.transAxes, fontsize=6.5, va="top",
                family="monospace",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#CCCCCC", alpha=0.9))

        ax.set_xlabel("PSI", fontsize=8)
        ax.set_xlim(0, 1)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.grid(True, lw=0.4, alpha=0.4, zorder=0)
        ax.set_axisbelow(True)

    axes[0].set_ylabel("Fraction of labeled bins", fontsize=8)

    # Shared legend for PSI categories
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(color="#5E9ACC", label="PSI < 0.1 (skipped)"),
        Patch(color="#E67E22", label="PSI 0.1–0.9 (regulated)"),
        Patch(color="#C0392B", label="PSI > 0.9 (included)"),
    ]
    fig.legend(handles=legend_handles, fontsize=7, loc="upper right",
               bbox_to_anchor=(1.0, 1.0), framealpha=0.9, edgecolor="#CCCCCC")

    fig.suptitle("Genome-wide PSI label distribution (masked bins only)", fontsize=9, y=1.02)
    plt.tight_layout()

    out = os.path.join(OUT_DIR, "panel_E_distribution.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    print(f"Saved: {out}")

    print("\n=== PSI category breakdown ===")
    for ct in CELL_TYPES:
        vals = psi_per_ct[ct]
        lo, mid, hi = category_fractions(vals)
        print(f"  {ct}: {len(vals):,} masked bins | "
              f"<0.1: {lo:.1f}% | 0.1-0.9: {mid:.1f}% | >0.9: {hi:.1f}%")


if __name__ == "__main__":
    print(f"Loading up to {MAX_SHARDS} train shards from:\n  {SHARD_DIR}\n")
    psi_per_ct = load_psi_values(SHARD_DIR)
    plot(psi_per_ct)
