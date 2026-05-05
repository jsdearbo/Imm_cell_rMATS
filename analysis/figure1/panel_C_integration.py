"""
panel_C_integration.py

Figure 1C: rMATS PSI integration — showing how the binary skeleton is
corrected with measured PSI values at splicing event loci.

Same Syk locus as Panel B. Three tracks stacked:
  1. Gene model (same as B)
  2. Binary skeleton (before rMATS integration)
  3. PSI-corrected label (after rMATS integration)

The skipped exon (52596957–52597054) is colored by PSI (0.22 in B cells)
in the bottom track, illustrating the overwrite.

Cell type: B Cell (PSI values from B_only SE.MATS.JC.txt)

Output: figure1_out/panel_C_integration.pdf
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

# ── Paths ─────────────────────────────────────────────────────────────────────
GTF_PATH = os.path.join(
    os.path.dirname(__file__),
    "../../results/stringtie_per_rep/dominant_isoforms/B_Cell_dominant.gtf",
)
SE_PATH = os.path.join(
    os.path.dirname(__file__),
    "../../results/B_only/SE.MATS.JC.txt",
)
OUT_DIR = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Locus (shared with Panel B) ───────────────────────────────────────────────
CHROM      = "chr13"
VIEW_START = 52_582_000
VIEW_END   = 52_613_000
GENE_NAME  = "Syk"
STRAND     = "+"
BIN_SIZE   = 32

# PSI colormap: blue (low PSI) → white → red (high PSI)
PSI_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "psi_bwr", ["#2166AC", "#F7F7F7", "#C0392B"]
)


def parse_gtf_exons(gtf_path, chrom, gene_name, view_start, view_end):
    exons = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[0] != chrom or parts[2] != "exon":
                continue
            if f'gene_name "{gene_name}"' not in parts[8]:
                continue
            s = max(int(parts[3]), view_start)
            e = min(int(parts[4]), view_end)
            if s < e:
                exons.append((s, e))
    return sorted(exons)


def make_binary_track(exons, view_start, view_end, bin_size):
    n_bins = (view_end - view_start) // bin_size
    track  = np.zeros(n_bins, dtype=float)
    for es, ee in exons:
        b0 = max(0, (es - view_start) // bin_size)
        b1 = min(n_bins, (ee - view_start + bin_size - 1) // bin_size)
        track[b0:b1] = 1.0
    return track


def load_se_events_in_view(se_path, chrom, view_start, view_end, cell_type_psi_col="IncLevel1"):
    """Return list of (exon_start, exon_end, mean_psi) for SE events in view."""
    df = pd.read_csv(se_path, sep="\t")
    df = df[df["chr"] == chrom].copy()
    df = df[(df["exonStart_0base"] >= view_start) & (df["exonEnd"] <= view_end)]

    def parse_psi(s):
        if pd.isna(s) or s == "NA":
            return np.nan
        vals = [float(x) for x in str(s).split(",") if x != "NA"]
        return float(np.mean(vals)) if vals else np.nan

    df["psi"] = df[cell_type_psi_col].apply(parse_psi)
    # Deduplicate by exon coordinates — keep highest-coverage event per exon
    df = df.dropna(subset=["psi"])
    df["ijc_total"] = df["IJC_SAMPLE_1"].apply(
        lambda s: sum(float(x) for x in str(s).split(",") if x not in ("", "NA"))
        if pd.notna(s) else 0.0
    )
    df = df.sort_values("ijc_total", ascending=False).drop_duplicates(
        subset=["exonStart_0base", "exonEnd"]
    )
    events = [(int(r["exonStart_0base"]), int(r["exonEnd"]), float(r["psi"]))
              for _, r in df.iterrows()]
    return events


def apply_psi_to_track(binary_track, se_events, view_start, bin_size):
    """Overwrite binary track at SE event coordinates with measured PSI values."""
    psi_track = binary_track.copy()
    n_bins = len(binary_track)
    for es, ee, psi in se_events:
        b0 = max(0, (es - view_start) // bin_size)
        b1 = min(n_bins, (ee - view_start + bin_size - 1) // bin_size)
        psi_track[b0:b1] = psi
    return psi_track


def draw_gene_model(ax, exons, view_start, view_end, strand, skip_exon=None):
    """Draw gene model track: backbone + exon boxes + optional highlighted exon."""
    x_lims = (view_start / 1e3, view_end / 1e3)
    ax.set_xlim(*x_lims)
    ax.set_ylim(-0.5, 1.5)
    ax.plot(x_lims, [0.5, 0.5], color="#555555", lw=1.0, zorder=1)

    n_ticks = max(1, (view_end - view_start) // 2000)
    for tp in np.linspace(view_start + 1000, view_end - 1000, n_ticks):
        dx = 500 if strand == "+" else -500
        ax.annotate("", xy=((tp + dx) / 1e3, 0.5), xytext=(tp / 1e3, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="#AAAAAA",
                                   lw=0.5, mutation_scale=6))
    for es, ee in exons:
        color = "#D6604D" if (skip_exon and es == skip_exon[0]) else "#2166AC"
        rect = mpatches.FancyBboxPatch(
            (es / 1e3, 0.1), (ee - es) / 1e3, 0.8,
            boxstyle="square,pad=0",
            facecolor=color, edgecolor="#1A4E7A" if color == "#2166AC" else "#A93226",
            lw=0.8, zorder=2,
        )
        ax.add_patch(rect)
    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)


def plot(exons, binary_track, psi_track, se_events, view_start, view_end, bin_size):
    skip_exon = se_events[0][:2] if se_events else None

    fig, axes = plt.subplots(
        3, 1, figsize=(7, 4.2),
        gridspec_kw={"height_ratios": [1.0, 0.75, 0.75]},
    )
    fig.subplots_adjust(hspace=0.06)

    n_bins = len(binary_track)
    x_kb   = view_start / 1e3 + (np.arange(n_bins) + 0.5) * bin_size / 1e3
    x_lims = (view_start / 1e3, view_end / 1e3)

    # ── Track 1: Gene model ───────────────────────────────────────────────────
    draw_gene_model(axes[0], exons, view_start, view_end, STRAND,
                    skip_exon=skip_exon)
    axes[0].set_title(
        f"{GENE_NAME} — rMATS PSI integration  ({CHROM})", fontsize=9
    )
    axes[0].set_ylabel("Gene model", fontsize=7.5, rotation=0,
                        labelpad=52, va="center")

    # ── Track 2: Binary skeleton ──────────────────────────────────────────────
    ax = axes[1]
    for xv, v in zip(x_kb, binary_track):
        ax.bar(xv, v, width=bin_size / 1e3,
               color="#2166AC" if v > 0.5 else "#CCCCCC", lw=0, zorder=2)
    ax.set_xlim(*x_lims)
    ax.set_ylim(-0.05, 1.2)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"], fontsize=7)
    ax.set_ylabel("Before\n(binary)", fontsize=7.5, rotation=0,
                   labelpad=52, va="center")
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False, labelsize=7)

    # ── Track 3: PSI-corrected label ──────────────────────────────────────────
    ax = axes[2]
    for xv, v in zip(x_kb, psi_track):
        if v < 0.01:
            c = "#CCCCCC"
        else:
            c = PSI_CMAP(v)
        ax.bar(xv, v, width=bin_size / 1e3, color=c, lw=0, zorder=2)

    ax.set_xlim(*x_lims)
    ax.set_ylim(-0.05, 1.2)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "0.5", "1"], fontsize=7)
    ax.set_ylabel("After\n(PSI label)", fontsize=7.5, rotation=0,
                   labelpad=52, va="center")
    ax.set_xlabel(f"Genomic position — {CHROM} (kb)", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)

    # Colorbar for PSI
    sm = plt.cm.ScalarMappable(cmap=PSI_CMAP, norm=mcolors.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[2], orientation="vertical",
                        fraction=0.015, pad=0.01, shrink=0.8)
    cbar.set_label("PSI", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # Annotate the skipped exon on track 3
    if skip_exon:
        es, ee = skip_exon
        mid_x = (es + ee) / 2 / 1e3
        psi_val = se_events[0][2]
        axes[2].annotate(
            f"PSI = {psi_val:.2f}",
            xy=(mid_x, psi_val + 0.05), xytext=(mid_x - 3.5, 0.85),
            fontsize=6.5, color="#C0392B", fontweight="bold",
            arrowprops=dict(arrowstyle="-|>", color="#C0392B", lw=0.7,
                            mutation_scale=7),
        )

    plt.tight_layout()
    out = os.path.join(OUT_DIR, "panel_C_integration.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    print(f"Saved: {out}")
    print(f"SE events in view: {len(se_events)}")
    for es, ee, psi in se_events:
        print(f"  exon {es}–{ee}: PSI = {psi:.3f}")


if __name__ == "__main__":
    exons     = parse_gtf_exons(GTF_PATH, CHROM, GENE_NAME, VIEW_START, VIEW_END)
    binary    = make_binary_track(exons, VIEW_START, VIEW_END, BIN_SIZE)
    se_events = load_se_events_in_view(SE_PATH, CHROM, VIEW_START, VIEW_END)
    psi_track = apply_psi_to_track(binary, se_events, VIEW_START, BIN_SIZE)
    plot(exons, binary, psi_track, se_events, VIEW_START, VIEW_END, BIN_SIZE)
