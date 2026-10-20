"""
panel_D_locus_example.py

Figure 1D: Example genomic locus showing PSI-labeled training window.
Scans train shards for a high-signal window with visible inter-cell-type
PSI variation, then renders a genome-browser-style figure.

Tracks (top → bottom):
  1. GC content (32-bp bins)
  2–4. Per-cell-type PSI (B Cell, Macrophage, T Cell) with confidence shading

Output:
  figure1_out/panel_D_locus_example.pdf   — full figure
  figure1_out/panel_D_locus_info.txt      — selected window coords + stats
"""

import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── Paths ─────────────────────────────────────────────────────────────────────
SHARD_DIR = os.path.join(
    os.path.dirname(__file__),
    "../../data_labeling/results/labelled_data/multitask_balanced/train",
)
OUT_DIR = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

CELL_TYPES  = ["B Cell", "Macrophage", "T Cell"]   # task indices 0, 1, 2
TASK_COLORS = {
    "B Cell":     "#2166AC",
    "Macrophage": "#4DAC26",
    "T Cell":     "#D6604D",
}

BIN_SIZE    = 32     # bp per bin
MAX_SHARDS  = 4      # shards to scan for best example (--scan mode only)

# Pre-selected best window (from scan over 8 shards).
# Pass --scan to rediscover; otherwise this is loaded directly.
DEFAULT_SHARD = "shard_10000.npz"
DEFAULT_IDX   = 70


def gc_content(one_hot, n_bins):
    """Bin sequence one-hot (L,4) → GC fraction per bin."""
    L = one_hot.shape[0]
    bp_per_bin = L // n_bins
    gc = (one_hot[:, 1] + one_hot[:, 2])          # C=1, G=2
    return gc[: n_bins * bp_per_bin].reshape(n_bins, bp_per_bin).mean(axis=1)


def score_window(targets, mask, confidence):
    """
    Score a window for visual interest:
    - Must have at least 50 masked bins in each cell type
    - Score = inter-cell-type std across masked bins (higher = more interesting)
    """
    for ti in range(3):
        if mask[:, ti].sum() < 50:
            return -1.0
    masked_any = mask.any(axis=1)
    psi_where_masked = targets[masked_any]
    return float(psi_where_masked.std(axis=0).mean())


def find_best_example(shard_dir, max_shards=MAX_SHARDS):
    """
    Two-pass approach:
    1. Scan small arrays (targets, mask, confidence, has_signal, coords) with mmap
       to find the best-scoring window without loading sequence into RAM.
    2. Load the full shard for the winner to get the sequence array.
    """
    shards = sorted(glob.glob(os.path.join(shard_dir, "*.npz")))[:max_shards]
    best_score  = -1.0
    best_idx    = None
    best_shard  = None

    SMALL_KEYS = {"targets", "mask", "confidence", "has_signal",
                  "chrom", "window_start", "window_end"}

    for path in shards:
        d = np.load(path, allow_pickle=True, mmap_mode="r")
        n_examples = d["targets"].shape[0]
        for i in range(n_examples):
            if not d["has_signal"][i]:
                continue
            score = score_window(d["targets"][i], d["mask"][i], d["confidence"][i])
            if score > best_score:
                best_score = score
                best_idx   = i
                best_shard = path
        print(f"  scanned {os.path.basename(path)}")

    # Full load for the winning shard to get sequence
    d = np.load(best_shard, allow_pickle=True)
    best_data = {
        "targets":      d["targets"][best_idx],
        "mask":         d["mask"][best_idx],
        "confidence":   d["confidence"][best_idx],
        "sequence":     d["sequence"][best_idx],
        "chrom":        str(d["chrom"][best_idx]),
        "window_start": int(d["window_start"][best_idx]),
        "window_end":   int(d["window_end"][best_idx]),
    }
    print(f"\nBest window: {best_data['chrom']}:{best_data['window_start']}-{best_data['window_end']}"
          f"  (score={best_score:.4f}, shard={os.path.basename(best_shard)}, idx={best_idx})")
    return best_data


def plot(data):
    targets    = data["targets"]     # (bins, 3)
    mask       = data["mask"]        # (bins, 3)
    confidence = data["confidence"]  # (bins, 3)
    seq        = data["sequence"]    # (L, 4)
    chrom      = data["chrom"]
    w_start    = data["window_start"]
    w_end      = data["window_end"]

    n_bins = targets.shape[0]
    gc     = gc_content(seq, n_bins)

    # Genomic x-axis in Mb
    bp_centers = w_start + (np.arange(n_bins) + 0.5) * BIN_SIZE
    x_mb       = bp_centers / 1e6

    n_tracks = 1 + len(CELL_TYPES)
    fig, axes = plt.subplots(
        n_tracks, 1,
        figsize=(9, 1.7 * n_tracks),
        sharex=True,
        gridspec_kw={"height_ratios": [0.6] + [1.0] * len(CELL_TYPES)},
    )
    fig.subplots_adjust(hspace=0.05)

    # ── GC track ──────────────────────────────────────────────────────────────
    ax = axes[0]
    ax.fill_between(x_mb, gc, alpha=0.55, color="#888888", lw=0)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "0.5", "1"], fontsize=6)
    ax.set_ylabel("GC", fontsize=7, rotation=0, labelpad=28, va="center")
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.tick_params(bottom=False, labelsize=6)

    # ── PSI tracks ────────────────────────────────────────────────────────────
    for ti, ct in enumerate(CELL_TYPES):
        ax = axes[ti + 1]
        color = TASK_COLORS[ct]

        psi  = targets[:, ti].copy()
        msk  = mask[:, ti].astype(bool)
        conf = confidence[:, ti].copy()

        # Unmasked regions drawn as flat gray
        ax.fill_between(x_mb, np.where(~msk, 0, np.nan),
                        np.where(~msk, 1, np.nan),
                        color="#F0F0F0", lw=0, zorder=0)

        # PSI signal where masked
        psi_masked = np.where(msk, psi, np.nan)
        ax.fill_between(x_mb, psi_masked, alpha=0.35, color=color, lw=0, zorder=1)
        ax.plot(x_mb, psi_masked, color=color, lw=0.8, zorder=2)

        ax.set_ylim(-0.05, 1.05)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "", "1"], fontsize=6)
        ax.set_ylabel(ct, fontsize=7, rotation=0, labelpad=28, va="center",
                      color=color, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=6)
        if ti < len(CELL_TYPES) - 1:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(bottom=False)

    # ── x-axis label ──────────────────────────────────────────────────────────
    axes[-1].set_xlabel(f"Genomic position — {chrom} (Mb)", fontsize=8)
    axes[-1].xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    axes[-1].tick_params(axis="x", labelsize=7)

    fig.suptitle(
        f"PSI label tracks — {chrom}:{w_start:,}–{w_end:,}",
        fontsize=9, y=0.99,
    )

    out = os.path.join(OUT_DIR, "panel_D_locus_example.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    print(f"Saved: {out}")

    info_out = os.path.join(OUT_DIR, "panel_D_locus_info.txt")
    with open(info_out, "w") as fh:
        fh.write(f"chrom={chrom}\n")
        fh.write(f"window_start={w_start}\n")
        fh.write(f"window_end={w_end}\n")
        for ti, ct in enumerate(CELL_TYPES):
            msk = mask[:, ti].astype(bool)
            psi_m = targets[:, ti][msk]
            fh.write(f"{ct}: {msk.sum()} masked bins, "
                     f"mean PSI={psi_m.mean():.3f}, std={psi_m.std():.3f}\n")
    print(f"Saved: {info_out}")


def load_specific(shard_dir, shard_name, idx):
    """Load a specific (shard, index) pair directly — fast path."""
    path = os.path.join(shard_dir, shard_name)
    d    = np.load(path, allow_pickle=True)
    return {
        "targets":      d["targets"][idx],
        "mask":         d["mask"][idx],
        "confidence":   d["confidence"][idx],
        "sequence":     d["sequence"][idx],
        "chrom":        str(d["chrom"][idx]),
        "window_start": int(d["window_start"][idx]),
        "window_end":   int(d["window_end"][idx]),
    }


if __name__ == "__main__":
    import sys
    if "--scan" in sys.argv:
        print(f"Scanning up to {MAX_SHARDS} shards for best example...")
        data = find_best_example(SHARD_DIR)
    else:
        print(f"Loading pre-selected example: {DEFAULT_SHARD} idx={DEFAULT_IDX}")
        print("  (pass --scan to rediscover the best window)")
        data = load_specific(SHARD_DIR, DEFAULT_SHARD, DEFAULT_IDX)
        print(f"  window: {data['chrom']}:{data['window_start']}-{data['window_end']}")
    plot(data)
