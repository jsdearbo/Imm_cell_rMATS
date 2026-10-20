"""
verify_npz_labels.py

Sanity-check: confirm that PSI values shown in panel C (from master BEDs)
are consistent with what is actually stored in the NPZ training shards.

Three sources are compared at each locus:
  1. Master BED   — what panel C plots directly
  2. Per-sample avg — replicates the NPZ-build averaging (MIN_READS filter)
  3. NPZ extracted  — the actual values the model trained on

Output: figure1_out/verify_npz_{ptprc,pkm}.pdf  +  printed stats table
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = os.path.join(os.path.dirname(__file__), "../..")
MASTER_BED = {
    "B_Cell":     os.path.join(BASE, "data_labeling/results/PSI_labelled_bed/master_beds/B_Cell_master.bed"),
    "Macrophage": os.path.join(BASE, "data_labeling/results/PSI_labelled_bed/master_beds/Macrophage_master.bed"),
    "T_Cell":     os.path.join(BASE, "data_labeling/results/PSI_labelled_bed/master_beds/T_Cell_master.bed"),
}
MANIFEST_PATH = os.path.join(BASE, "data_labeling/results/labelled_data/multitask_balanced/manifest.tsv")
NPZ_ROOT      = os.path.join(BASE, "data_labeling/results/labelled_data/multitask_balanced")
OUT_DIR       = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────────────
BIN_SIZE = 32
MIN_READS = 10  # matches data_conversion_fast.py:MIN_READS_PER_REPLICATE
NPZ_TASK_ORDER = ["B_Cell", "Macrophage", "T_Cell"]  # index 0, 1, 2

LOCI = {
    "ptprc": {
        "chrom":      "chr1",
        "view_start": 138_113_300,
        "view_end":   138_126_526,
        "highlight":  [(138_122_456, 138_122_584),
                       (138_120_020, 138_120_166),
                       (138_117_754, 138_117_894)],
    },
    "pkm": {
        "chrom":      "chr9",
        "view_start":  59_673_000,
        "view_end":    59_680_500,
        "highlight":  [(59_675_028, 59_675_195)],
    },
}

CELL_COLORS = {"B_Cell": "#2166AC", "Macrophage": "#4DAC26", "T_Cell": "#D6604D"}

# ── Rasterization ─────────────────────────────────────────────────────────────

def rasterize_bed(bed_rows, view_start, view_end, bin_size=BIN_SIZE):
    """
    Length-weighted PSI per bin, matching data_conversion_fast.py logic.
    bed_rows: list of (start, end, psi, conf)
    Returns: (psi_bins, conf_bins, event_bins) each of shape (n_bins,)
    """
    n_bins = (view_end - view_start) // bin_size
    psi_acc   = np.zeros(n_bins, dtype=np.float64)
    conf_acc  = np.zeros(n_bins, dtype=np.float64)
    event_acc = np.zeros(n_bins, dtype=np.float64)

    for s, e, psi, conf in bed_rows:
        if np.isnan(psi):
            continue
        s_in = max(s, view_start)
        e_in = min(e, view_end)
        if s_in >= e_in:
            continue
        b0 = (s_in - view_start) // bin_size
        b1 = min((e_in - view_start + bin_size - 1) // bin_size, n_bins)
        for b in range(int(b0), int(b1)):
            bin_s = view_start + b * bin_size
            bin_e = bin_s + bin_size
            ov = max(0, min(e_in, bin_e) - max(s_in, bin_s))
            if ov <= 0:
                continue
            psi_acc[b]   += ov * psi
            conf_acc[b]  += ov * conf
            event_acc[b] += ov

    nz = event_acc > 0
    out_psi  = np.where(nz, psi_acc  / event_acc, np.nan)
    out_conf = np.where(nz, conf_acc / event_acc, np.nan)
    return out_psi, out_conf, event_acc


def load_bed(path, chrom, view_start, view_end):
    """Load a 5-col BED (no header: chrom start end psi conf), filter to region."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chrom", "start", "end", "psi", "conf"],
                     dtype={"chrom": str, "start": int, "end": int,
                            "psi": float, "conf": float})
    df = df[(df["chrom"] == chrom) &
            (df["end"] > view_start) & (df["start"] < view_end)]
    return list(zip(df["start"], df["end"], df["psi"], df["conf"]))


# ── Master BED rasterization ──────────────────────────────────────────────────

def master_bed_raster(cell_type, chrom, view_start, view_end):
    rows = load_bed(MASTER_BED[cell_type], chrom, view_start, view_end)
    psi, conf, events = rasterize_bed(rows, view_start, view_end)
    return psi, events


# ── Per-sample average (mirrors data_conversion_fast.py) ─────────────────────

def per_sample_avg(manifest_df, cell_type, chrom, view_start, view_end):
    """
    Average per-replicate rasters, applying MIN_READS_PER_REPLICATE filter
    on rMATS-derived bins (skeleton bins: conf==0.0, always included).
    """
    n_bins = (view_end - view_start) // BIN_SIZE
    acc_targ = np.zeros(n_bins, dtype=np.float64)
    acc_mask = np.zeros(n_bins, dtype=np.float64)

    rep_paths = manifest_df[manifest_df["cell_type"] == cell_type]["bed_path"].tolist()
    for path in rep_paths:
        rows = load_bed(path, chrom, view_start, view_end)
        psi_r, conf_r, events_r = rasterize_bed(rows, view_start, view_end)

        if not np.any(events_r > 0):
            continue

        reads = np.expm1(conf_r)
        is_skeleton = (conf_r == 0.0)
        rep_mask = ((events_r > 0) & (is_skeleton | (reads >= MIN_READS))).astype(float)

        psi_r_clean = np.where(np.isnan(psi_r), 0.0, psi_r)
        acc_targ += psi_r_clean * rep_mask
        acc_mask += rep_mask

    with np.errstate(divide="ignore", invalid="ignore"):
        avg = np.where(acc_mask > 0, acc_targ / acc_mask, np.nan)
    return avg, acc_mask


# ── NPZ extraction ────────────────────────────────────────────────────────────

def extract_from_npz(chrom, view_start, view_end):
    """
    Scan all shards in train/val/test for windows fully containing the locus.
    Returns dict keyed by cell_type → (psi_bins, mask_bins) at 32-bp resolution.
    If multiple windows overlap (different strides), average masked PSI across them.
    """
    n_bins = (view_end - view_start) // BIN_SIZE
    acc = {ct: {"psi": np.zeros(n_bins), "mask": np.zeros(n_bins), "count": np.zeros(n_bins)}
           for ct in NPZ_TASK_ORDER}
    hits = 0

    for split in ("train", "val", "test"):
        shard_dir = os.path.join(NPZ_ROOT, split)
        if not os.path.isdir(shard_dir):
            continue
        for shard_path in sorted(glob.glob(os.path.join(shard_dir, "shard_*.npz"))):
            data = np.load(shard_path, allow_pickle=False)
            chroms       = data["chrom"]
            win_starts   = data["window_start"]
            win_ends     = data["window_end"]

            for i in range(len(chroms)):
                if chroms[i] != chrom:
                    continue
                ws, we = int(win_starts[i]), int(win_ends[i])
                # Window must fully contain the locus
                if ws > view_start or we < view_end:
                    continue
                hits += 1

                # Bin offsets within the window (force exact n_bins length)
                b0 = (view_start - ws) // BIN_SIZE
                b1 = b0 + n_bins

                targets    = data["targets"][i]    # (n_window_bins, 3)
                masks      = data["mask"][i]       # (n_window_bins, 3)

                for ti, ct in enumerate(NPZ_TASK_ORDER):
                    psi_slice  = targets[b0:b1, ti]
                    mask_slice = masks[b0:b1, ti]
                    acc[ct]["psi"]   += psi_slice * mask_slice
                    acc[ct]["mask"]  += mask_slice
                    acc[ct]["count"] += (mask_slice > 0).astype(float)

            data.close()

    if hits == 0:
        print(f"  [WARN] No NPZ windows found fully containing {chrom}:{view_start}-{view_end}")

    result = {}
    for ct in NPZ_TASK_ORDER:
        m = acc[ct]["mask"]
        psi = np.where(m > 0, acc[ct]["psi"] / m, np.nan)
        result[ct] = (psi, (m > 0).astype(float))

    print(f"  NPZ: {hits} overlapping window(s) found for {chrom}:{view_start}-{view_end}")
    return result


# ── Pearson r helper ──────────────────────────────────────────────────────────

def pearsonr_masked(a, b, mask):
    valid = mask > 0
    if valid.sum() < 3:
        return float("nan"), float("nan")
    av, bv = a[valid], b[valid]
    if av.std() < 1e-8 or bv.std() < 1e-8:
        return float("nan"), float("nan")
    r = np.corrcoef(av, bv)[0, 1]
    mae = np.abs(av - bv).mean()
    return float(r), float(mae)


# ── Plot one locus ────────────────────────────────────────────────────────────

def plot_locus(locus_key):
    locus      = LOCI[locus_key]
    chrom      = locus["chrom"]
    view_start = locus["view_start"]
    view_end   = locus["view_end"]
    hi_exons   = locus.get("highlight", [])
    n_bins     = (view_end - view_start) // BIN_SIZE
    x_kb       = view_start / 1e3 + (np.arange(n_bins) + 0.5) * BIN_SIZE / 1e3
    x_lims     = (view_start / 1e3, view_end / 1e3)

    print(f"\n{'='*60}")
    print(f"Locus: {locus_key.upper()}  ({chrom}:{view_start:,}–{view_end:,})")
    print(f"{'='*60}")

    manifest_df = pd.read_csv(MANIFEST_PATH, sep="\t", dtype=str)

    print("Loading master BEDs...")
    master = {}
    for ct in NPZ_TASK_ORDER:
        psi, events = master_bed_raster(ct, chrom, view_start, view_end)
        master[ct] = (psi, events)
        n_valid = int(np.sum(~np.isnan(psi)))
        print(f"  {ct}: {n_valid} non-NaN bins")

    print("Computing per-sample averages...")
    persample = {}
    for ct in NPZ_TASK_ORDER:
        avg, mask = per_sample_avg(manifest_df, ct, chrom, view_start, view_end)
        persample[ct] = (avg, mask)
        n_valid = int(np.sum(~np.isnan(avg)))
        print(f"  {ct}: {n_valid} non-NaN bins")

    print("Scanning NPZ shards...")
    npz = extract_from_npz(chrom, view_start, view_end)

    # ── Print stats ────────────────────────────────────────────────────────────
    print(f"\n{'─'*60}")
    print(f"{'Cell type':<14} {'Source pair':<28} {'Pearson r':>10} {'MAE':>8}")
    print(f"{'─'*60}")
    for ct in NPZ_TASK_ORDER:
        npz_psi, npz_mask = npz[ct]
        npz_psi_clean = np.where(np.isnan(npz_psi), 0.0, npz_psi)

        # master BED vs NPZ
        m_psi, _ = master[ct]
        m_psi_clean = np.where(np.isnan(m_psi), 0.0, m_psi)
        r, mae = pearsonr_masked(m_psi_clean, npz_psi_clean, npz_mask)
        print(f"{ct:<14} {'Master BED vs NPZ':<28} {r:>10.4f} {mae:>8.4f}")

        # per-sample avg vs NPZ
        ps_psi, _ = persample[ct]
        ps_psi_clean = np.where(np.isnan(ps_psi), 0.0, ps_psi)
        r2, mae2 = pearsonr_masked(ps_psi_clean, npz_psi_clean, npz_mask)
        print(f"{ct:<14} {'Per-sample avg vs NPZ':<28} {r2:>10.4f} {mae2:>8.4f}")
    print(f"{'─'*60}")

    # ── Figure: 3 cell types × 3 sources ──────────────────────────────────────
    sources = ["Master BED\n(panel C)", "Per-sample avg\n(NPZ source)", "NPZ extracted\n(model trained on)"]
    fig, axes = plt.subplots(
        3, 3, figsize=(13, 7), sharex=True, sharey=True,
        gridspec_kw={"hspace": 0.05, "wspace": 0.05},
    )
    fig.suptitle(f"{locus_key.upper()}  ({chrom}:{view_start:,}–{view_end:,})  —  PSI label consistency",
                 fontsize=10, y=0.98)

    for col_i, ct in enumerate(NPZ_TASK_ORDER):
        color = CELL_COLORS[ct]
        npz_psi, npz_mask = npz[ct]

        tracks = [
            (master[ct][0],     master[ct][1] > 0),   # row 0
            (persample[ct][0],  persample[ct][1] > 0), # row 1
            (npz_psi,           npz_mask > 0),          # row 2
        ]

        for row_i, (psi, mask) in enumerate(tracks):
            ax = axes[row_i][col_i]
            # Background
            ax.fill_between(x_lims, [0, 0], [0.02, 0.02], color="#EEEEEE", lw=0)
            # Bars — gray where mask==0, colored where mask==1
            bar_colors = np.where(mask, color, "#CCCCCC")
            psi_plot = np.where(np.isnan(psi), 0.0, psi)
            ax.bar(x_kb, np.maximum(psi_plot, 0.01), width=BIN_SIZE / 1e3,
                   color=bar_colors, lw=0, zorder=2)

            ax.set_xlim(*x_lims)
            ax.set_ylim(-0.05, 1.20)
            ax.set_yticks([0, 0.5, 1.0])
            ax.tick_params(labelsize=6)
            ax.spines[["top", "right"]].set_visible(False)

            # Highlight exons
            for hs, he in hi_exons:
                ax.axvspan(hs / 1e3, he / 1e3, color="#FFE0B2", alpha=0.35, zorder=0)

            # Row labels (left column only)
            if col_i == 0:
                ax.set_ylabel(sources[row_i], fontsize=7, rotation=0,
                              labelpad=80, va="center")

            # Column headers (top row only)
            if row_i == 0:
                ax.set_title(ct.replace("_", " "), fontsize=8, color=color, pad=4)

            # x-axis label (bottom row only)
            if row_i == 2:
                ax.set_xlabel(f"{chrom} (kb)", fontsize=7)
                ax.tick_params(axis="x", labelsize=6)
            else:
                ax.tick_params(labelbottom=False)

            # Stats annotation
            if row_i < 2:
                npz_psi_clean = np.where(np.isnan(npz_psi), 0.0, npz_psi)
                psi_clean = np.where(np.isnan(psi), 0.0, psi)
                r, mae = pearsonr_masked(psi_clean, npz_psi_clean, npz_mask)
                if not np.isnan(r):
                    ax.text(0.98, 0.96, f"r={r:.3f}\nMAE={mae:.3f}",
                            transform=ax.transAxes, fontsize=5.5,
                            ha="right", va="top", color="#333333",
                            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))

    out = os.path.join(OUT_DIR, f"verify_npz_{locus_key}.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"\nSaved: {out}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    targets = sys.argv[1:] if len(sys.argv) > 1 else list(LOCI.keys())
    for key in targets:
        if key not in LOCI:
            print(f"Unknown locus '{key}'. Available: {list(LOCI.keys())}")
            continue
        plot_locus(key)
