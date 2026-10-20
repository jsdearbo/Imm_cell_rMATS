"""
panel_B_skeleton.py

Figure 1B: StringTie dominant-isoform binary label skeleton.
Shows how the dominant transcript annotation provides a genome-wide
scaffold where exon bins = 1 and intron/intergenic bins = 0.

Locus: Syk, chr13:52582000–52612000 (B Cell dominant isoform)
  — captures the first 3 exons and the alternatively-skipped exon (exon 2b)
  — the skipped exon (52596957–52597054, PSI=0.22 in B cells) appears
    here with label = 1 (skeleton), showing it before PSI correction

Output: figure1_out/panel_B_skeleton.pdf
"""

import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ─────────────────────────────────────────────────────────────────────
GTF_PATH = os.path.join(
    os.path.dirname(__file__),
    "../../results/stringtie_per_rep/dominant_isoforms/B_Cell_dominant.gtf",
)
OUT_DIR = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Locus ─────────────────────────────────────────────────────────────────────
CHROM      = "chr13"
VIEW_START = 52_582_000
VIEW_END   = 52_613_000
GENE_NAME  = "Syk"
STRAND     = "+"
BIN_SIZE   = 32   # bp per label bin


def parse_gtf_exons(gtf_path, chrom, gene_name, view_start, view_end):
    """Return list of (exon_start, exon_end) for gene_name on chrom within view."""
    exons = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[0] != chrom or parts[2] != "exon":
                continue
            attrs = parts[8]
            if f'gene_name "{gene_name}"' not in attrs:
                continue
            s, e = int(parts[3]), int(parts[4])
            # Restrict to view region
            s = max(s, view_start)
            e = min(e, view_end)
            if s < e:
                exons.append((s, e))
    return sorted(exons)


def make_binary_track(exons, view_start, view_end, bin_size):
    """Return 0/1 array at bin_size resolution covering [view_start, view_end)."""
    n_bins = (view_end - view_start) // bin_size
    track = np.zeros(n_bins, dtype=float)
    for es, ee in exons:
        b_start = max(0, (es - view_start) // bin_size)
        b_end   = min(n_bins, (ee - view_start + bin_size - 1) // bin_size)
        track[b_start:b_end] = 1.0
    return track


def plot(exons, binary_track, view_start, view_end, bin_size):
    fig, axes = plt.subplots(
        2, 1, figsize=(7, 3.0),
        gridspec_kw={"height_ratios": [1.2, 0.8]},
    )
    fig.subplots_adjust(hspace=0.08)

    n_bins = len(binary_track)
    x_kb   = view_start / 1e3 + (np.arange(n_bins) + 0.5) * bin_size / 1e3
    x_lims = (view_start / 1e3, view_end / 1e3)

    # ── Track 1: Gene model ───────────────────────────────────────────────────
    ax = axes[0]
    ax.set_xlim(*x_lims)
    ax.set_ylim(-0.5, 1.5)

    # Intron backbone
    ax.plot(x_lims, [0.5, 0.5], color="#555555", lw=1.0, zorder=1)

    # Strand arrow ticks
    n_ticks = max(1, (view_end - view_start) // 2000)
    tick_positions = np.linspace(view_start + 1000, view_end - 1000, n_ticks)
    for tp in tick_positions:
        dx = 500 if STRAND == "+" else -500
        ax.annotate("", xy=((tp + dx) / 1e3, 0.5),
                    xytext=(tp / 1e3, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="#AAAAAA",
                                   lw=0.5, mutation_scale=6))

    # Exon boxes
    for es, ee in exons:
        width = (ee - es) / 1e3
        rect  = mpatches.FancyBboxPatch(
            (es / 1e3, 0.1), width, 0.8,
            boxstyle="square,pad=0",
            facecolor="#2166AC", edgecolor="#1A4E7A", lw=0.8, zorder=2
        )
        ax.add_patch(rect)

    ax.set_ylabel("Gene model", fontsize=7.5, rotation=0, labelpad=52, va="center")
    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)
    ax.set_title(f"{GENE_NAME} — binary label skeleton  ({CHROM})", fontsize=9)

    # Annotation for the skipped exon
    skip_es, skip_ee = 52_596_957, 52_597_054
    if VIEW_START <= skip_es <= VIEW_END:
        mid_x = (skip_es + skip_ee) / 2 / 1e3
        ax.annotate("skipped exon\n(PSI = 0.22 in B cells)",
                    xy=(mid_x, 0.95), xytext=(mid_x, 1.4),
                    fontsize=6, ha="center", color="#C0392B",
                    arrowprops=dict(arrowstyle="-", color="#C0392B", lw=0.7))

    # ── Track 2: Binary label ─────────────────────────────────────────────────
    ax = axes[1]

    exon_color    = "#2166AC"
    intron_color  = "#CCCCCC"

    colors_per_bin = np.where(binary_track > 0.5, exon_color, intron_color)
    for xi, (xv, c) in enumerate(zip(x_kb, colors_per_bin)):
        ax.bar(xv, binary_track[xi], width=bin_size / 1e3,
               color=c, lw=0, zorder=2)

    ax.set_xlim(*x_lims)
    ax.set_ylim(-0.05, 1.2)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"], fontsize=7)
    ax.set_ylabel("PSI label\n(binary skeleton)", fontsize=7.5,
                  rotation=0, labelpad=52, va="center")
    ax.set_xlabel(f"Genomic position — {CHROM} (kb)", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)

    # Legend
    from matplotlib.patches import Patch
    ax.legend(
        handles=[Patch(color=exon_color, label="Exon (label = 1)"),
                 Patch(color=intron_color, label="Intron (label = 0)")],
        fontsize=6.5, loc="upper right", framealpha=0.9, edgecolor="#CCCCCC",
    )

    plt.tight_layout()
    out = os.path.join(OUT_DIR, "panel_B_skeleton.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    print(f"Saved: {out}")
    print(f"Exons in view: {len(exons)}")
    print(f"Binary track: {int(binary_track.sum())} / {len(binary_track)} bins = 1")


if __name__ == "__main__":
    exons  = parse_gtf_exons(GTF_PATH, CHROM, GENE_NAME, VIEW_START, VIEW_END)
    track  = make_binary_track(exons, VIEW_START, VIEW_END, BIN_SIZE)
    plot(exons, track, VIEW_START, VIEW_END, BIN_SIZE)
