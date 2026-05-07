"""
panel_C_integration.py

Figure 1C: RNA-seq coverage + PSI training labels at representative alternatively-spliced loci.

Two example loci (separate PDFs), each showing 10 tracks on a shared genomic x-axis:
  Rows 0–2: merged RNA-seq coverage per cell type (bigwig, linear scale)
  Row  3:   gene model (exon boxes labeled E{n}, intron lines, strand arrows)
  Rows 4–6: constitutive skeleton — dominant isoform binary labels at 32-bp bins
  Rows 7–9: PSI-corrected training labels per cell type (master BED intervals)

PSI labels come directly from the master BED files — the source used to build the NPZ shards.

Output: figure1_out/panel_C_{ptprc,pkm}.pdf
"""

import os
import re
import sys
import numpy as np
import pandas as pd
import pyBigWig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE      = os.path.join(os.path.dirname(__file__), "../..")
GTF_FILES = {
    "B Cell":     os.path.join(BASE, "results/stringtie_per_rep/dominant_isoforms/B_Cell_dominant.gtf"),
    "Macrophage": os.path.join(BASE, "results/stringtie_per_rep/dominant_isoforms/Macrophage_dominant.gtf"),
    "T Cell":     os.path.join(BASE, "results/stringtie_per_rep/dominant_isoforms/T_Cell_dominant.gtf"),
}
BED_FILES = {
    "B Cell":     os.path.join(BASE, "data_labeling/results/PSI_labelled_bed/master_beds/B_Cell_master.bed"),
    "Macrophage": os.path.join(BASE, "data_labeling/results/PSI_labelled_bed/master_beds/Macrophage_master.bed"),
    "T Cell":     os.path.join(BASE, "data_labeling/results/PSI_labelled_bed/master_beds/T_Cell_master.bed"),
}
BW_FILES = {
    "B Cell":     os.path.join(BASE, "results/merged_bw/B_Cell_merged.bw"),
    "Macrophage": os.path.join(BASE, "results/merged_bw/Macrophage_merged.bw"),
    "T Cell":     os.path.join(BASE, "results/merged_bw/T_Cell_merged.bw"),
}
OUT_DIR = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

BIN_SIZE = 32  # bp per label bin

# ── Visual constants ───────────────────────────────────────────────────────────
CELL_TYPES = ["B Cell", "Macrophage", "T Cell"]
CELL_COLORS = {
    "B Cell":     "#2166AC",
    "Macrophage": "#4DAC26",
    "T Cell":     "#D6604D",
}

# ── Locus definitions ─────────────────────────────────────────────────────────
LOCI = {
    "ptprc": {
        "label":      "Ptprc",
        "gene_name":  "Ptprc",
        "chrom":      "chr1",
        "strand":     "-",
        # exons 3–9 region + flank (exon 3 at 138126326, exon 9 at 138113564)
        "view_start": 138_113_300,
        "view_end":   138_126_526,
        # Exons 4, 5, 6 are the CD45 variable exons (alt-spliced across cell types)
        "highlight_exons": [
            (138_122_456, 138_122_584),   # E4
            (138_120_020, 138_120_166),   # E5
            (138_117_754, 138_117_894),   # E6
        ],
        # Canonical CD45 exon labels (GENCODE vM23 ENSMUST00000183301 numbering)
        "exon_labels": {
            (138_126_303, 138_126_326): "E3",
            (138_122_456, 138_122_584): "E4",
            (138_120_020, 138_120_166): "E5",
            (138_117_754, 138_117_894): "E6",
            (138_116_954, 138_117_025): "E7",
            (138_115_619, 138_115_642): "E8",
            (138_113_564, 138_113_779): "E9",
        },
    },
    "pkm": {
        "label":      "Pkm",
        "gene_name":  "Pkm",
        "chrom":      "chr9",
        "strand":     "+",
        "view_start": 59_673_000,
        "view_end":   59_680_500,
        "highlight_exons": [
            (59_675_028, 59_675_195),
        ],
    },
}


# ── Data helpers ──────────────────────────────────────────────────────────────

def parse_gtf_exons_with_numbers(gtf_path, chrom, gene_name, view_start, view_end):
    """Return [(start, end, exon_number), ...] for exons within view, sorted by start."""
    exons = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[0] != chrom or parts[2] != "exon":
                continue
            if f'gene_name "{gene_name}"' not in parts[8]:
                continue
            m = re.search(r'exon_number "(\d+)"', parts[8])
            enum = int(m.group(1)) if m else 0
            s = max(int(parts[3]), view_start)
            e = min(int(parts[4]), view_end)
            if s < e:
                key = (s, e)
                if key not in exons:
                    exons[key] = enum
    return sorted([(s, e, n) for (s, e), n in exons.items()], key=lambda x: x[0])


def parse_union_gtf_exons(gtf_paths, chrom, gene_name, view_start, view_end):
    """Union of exons across multiple GTFs — ensures alt-spliced exons appear in gene model."""
    exons = {}
    for path in gtf_paths:
        for s, e, n in parse_gtf_exons_with_numbers(path, chrom, gene_name, view_start, view_end):
            if (s, e) not in exons:
                exons[(s, e)] = n
    return sorted([(s, e, n) for (s, e), n in exons.items()], key=lambda x: x[0])


def make_binary_track(exons_flat, view_start, view_end):
    """Return 0/1 float array at BIN_SIZE resolution covering [view_start, view_end)."""
    n_bins = (view_end - view_start) // BIN_SIZE
    track  = np.zeros(n_bins, dtype=float)
    for es, ee in exons_flat:
        b_start = max(0, (es - view_start) // BIN_SIZE)
        b_end   = min(n_bins, (ee - view_start + BIN_SIZE - 1) // BIN_SIZE)
        track[b_start:b_end] = 1.0
    return track


def load_bed_intervals(bed_path, chrom, view_start, view_end):
    """Load master BED intervals overlapping the view region. Returns [(start, end, psi), ...]."""
    df = pd.read_csv(
        bed_path, sep="\t", header=None,
        names=["chrom", "start", "end", "psi", "conf"],
        dtype={"chrom": str, "start": int, "end": int, "psi": float, "conf": float},
    )
    df = df[df["chrom"] == chrom]
    df = df[(df["end"] > view_start) & (df["start"] < view_end)]
    return [(int(r["start"]), int(r["end"]), float(r["psi"])) for _, r in df.iterrows()]


def load_bw_values(bw_path, chrom, start, end):
    """Query bigwig at 1-bp resolution (linear scale). Returns (x_kb, values)."""
    bw   = pyBigWig.open(bw_path)
    vals = np.array(bw.values(chrom, start, end), dtype=float)
    bw.close()
    vals = np.nan_to_num(vals, nan=0.0)
    x_kb = (start + np.arange(len(vals))) / 1e3
    return x_kb, vals


def compute_psi_annotations(bed_dict, highlight_exons):
    """
    For each highlight exon, compute length-weighted mean PSI per cell type
    from overlapping BED intervals.
    Returns {(hs, he): {ct: mean_psi}}.
    """
    annots = {}
    for hs, he in highlight_exons:
        ct_vals = {}
        for ct, intervals in bed_dict.items():
            overlapping = [(s, e, psi) for s, e, psi in intervals if e > hs and s < he]
            if overlapping:
                total = sum(min(e, he) - max(s, hs) for s, e, _ in overlapping)
                if total > 0:
                    wpsi = sum((min(e, he) - max(s, hs)) * psi for s, e, psi in overlapping)
                    ct_vals[ct] = wpsi / total
        if ct_vals:
            annots[(hs, he)] = ct_vals
    return annots


# ── Drawing helpers ───────────────────────────────────────────────────────────

def draw_rna_track(ax, x_kb, values, color):
    ax.fill_between(x_kb, values, alpha=0.75, color=color, lw=0, zorder=2)
    y_max = max(float(values.max()) * 1.1, 0.1)
    ax.set_ylim(0, y_max)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=6)
    ax.yaxis.set_major_locator(plt.MaxNLocator(2, integer=True))


def draw_skeleton_track(ax, exons_with_nums, view_start, view_end, color):
    """Constitutive binary skeleton at BIN_SIZE resolution. Exon bins = cell color, intron = gray."""
    exons_flat = [(s, e) for s, e, _ in exons_with_nums]
    track = make_binary_track(exons_flat, view_start, view_end)
    n_bins = len(track)
    x_kb   = view_start / 1e3 + (np.arange(n_bins) + 0.5) * BIN_SIZE / 1e3
    for xi, val in enumerate(track):
        c = color if val > 0.5 else "#EEEEEE"
        ax.bar(x_kb[xi], max(val, 0.015), width=BIN_SIZE / 1e3, color=c, lw=0, zorder=2)
    ax.set_ylim(-0.05, 1.2)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"], fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=6)


def draw_gene_model(ax, exons_with_nums, view_start, view_end, strand,
                    highlight_exons=(), exon_labels=None):
    """
    exon_labels: optional dict {(start, end): "label"} for custom exon labeling.
    Falls back to GTF exon_number when absent.
    """
    ax.set_ylim(-0.4, 1.5)
    x_lo, x_hi = view_start / 1e3, view_end / 1e3
    ax.plot([x_lo, x_hi], [0.5, 0.5], color="#555555", lw=0.9, zorder=1)

    # Strand direction arrows
    n_ticks = max(1, (view_end - view_start) // 2000)
    for tp in np.linspace(view_start + 1000, view_end - 1000, n_ticks):
        dx = 400 if strand == "+" else -400
        ax.annotate("", xy=((tp + dx) / 1e3, 0.5), xytext=(tp / 1e3, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="#AAAAAA",
                                   lw=0.5, mutation_scale=6))

    hi_set      = {(es, ee) for es, ee in highlight_exons}
    exon_coords = [(s, e) for s, e, _ in exons_with_nums]

    # Intron labels between consecutive exons
    for i in range(len(exon_coords) - 1):
        _, e1 = exon_coords[i]
        s2, _ = exon_coords[i + 1]
        if s2 - e1 > 500:
            mid = (e1 + s2) / 2 / 1e3
            ax.text(mid, 0.52, "intron", ha="center", va="bottom",
                    fontsize=4.5, color="#888888")

    # Exon boxes with labels
    for es, ee, enum in exons_with_nums:
        is_hi = (es, ee) in hi_set or any((hs <= es and ee <= he) for hs, he in hi_set)
        fc = "#E8820C" if is_hi else "#888888"
        ec = "#C06000" if is_hi else "#555555"
        rect = mpatches.FancyBboxPatch(
            (es / 1e3, 0.1), (ee - es) / 1e3, 0.8,
            boxstyle="square,pad=0", facecolor=fc, edgecolor=ec, lw=0.7, zorder=2,
        )
        ax.add_patch(rect)
        if exon_labels is not None:
            label = exon_labels.get((es, ee), "")
        else:
            label = f"E{enum}" if enum else ""
        if label and (ee - es) > 80:
            ax.text((es + ee) / 2 / 1e3, 0.5, label,
                    ha="center", va="center", fontsize=5,
                    color="white", fontweight="bold", zorder=3)

    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)


def draw_psi_track(ax, intervals, view_start, view_end, color="#2166AC"):
    """Draw PSI labels from master BED intervals."""
    ax.fill_between(
        [view_start / 1e3, view_end / 1e3], [0.02, 0.02], color="#EEEEEE", lw=0, zorder=0
    )
    for start, end, psi in intervals:
        s_kb = max(start, view_start) / 1e3
        e_kb = min(end,   view_end)   / 1e3
        if s_kb >= e_kb:
            continue
        cx  = (s_kb + e_kb) / 2
        w   = e_kb - s_kb
        if psi < 0.01:
            c, h = "#CCCCCC", 0.015
        else:
            c, h = color, psi
        ax.bar(cx, h, width=w, color=c, lw=0, zorder=2)

    ax.set_xlim(view_start / 1e3, view_end / 1e3)
    ax.set_ylim(-0.05, 1.15)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "", "1"], fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=6)


def _section_label(ax, text):
    """Add a bold section header inside the top-left of an axes."""
    ax.text(0.01, 0.93, text, transform=ax.transAxes,
            fontsize=6.5, fontweight="bold", color="#333333",
            ha="left", va="top", clip_on=True,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1))


# ── Main plot function ────────────────────────────────────────────────────────

def plot_locus(locus_key):
    locus      = LOCI[locus_key]
    chrom      = locus["chrom"]
    gene_name  = locus["gene_name"]
    strand     = locus["strand"]
    view_start = locus["view_start"]
    view_end   = locus["view_end"]
    hi_exons   = locus.get("highlight_exons", [])

    print(f"\n── {locus['label']} ({chrom}:{view_start:,}–{view_end:,}) ──")

    # ── Load data ─────────────────────────────────────────────────────────────
    # Gene model uses union of all cell-type dominant isoforms so alt-spliced exons appear
    exons = parse_union_gtf_exons(
        list(GTF_FILES.values()), chrom, gene_name, view_start, view_end
    )
    print(f"  GTF exons in view (gene model, union): {len(exons)}")

    skel_exons = {}
    for ct in CELL_TYPES:
        skel_exons[ct] = parse_gtf_exons_with_numbers(
            GTF_FILES[ct], chrom, gene_name, view_start, view_end
        )
        print(f"  {ct} skeleton exons: {len(skel_exons[ct])}")

    bed_dict = {}
    for ct in CELL_TYPES:
        intervals = load_bed_intervals(BED_FILES[ct], chrom, view_start, view_end)
        bed_dict[ct] = intervals
        psi_vals = [p for _, _, p in intervals]
        if psi_vals:
            print(f"  {ct}: {len(intervals)} BED intervals  PSI {min(psi_vals):.2f}–{max(psi_vals):.2f}")
        else:
            print(f"  {ct}: 0 BED intervals in view")

    bw_data = {}
    for ct in CELL_TYPES:
        print(f"  Loading bigwig for {ct}...")
        bw_data[ct] = load_bw_values(BW_FILES[ct], chrom, view_start, view_end)

    psi_annotations = compute_psi_annotations(bed_dict, hi_exons)
    for (hs, he), ct_psi in psi_annotations.items():
        print(f"  Highlight {hs:,}–{he:,}: " + "  ".join(f"{ct}={v:.2f}" for ct, v in ct_psi.items()))

    # ── Layout: 10 tracks ─────────────────────────────────────────────────────
    # 0: gene model   1–3: RNA-seq   4–6: skeleton   7–9: PSI labels
    height_ratios = [0.7, 0.8, 0.8, 0.8, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0]
    fig, axes = plt.subplots(
        10, 1,
        figsize=(7, sum(height_ratios) + 0.5),
        sharex=True,
        gridspec_kw={"height_ratios": height_ratios},
    )
    fig.subplots_adjust(hspace=0.25)

    # ── Gene model (row 0) ────────────────────────────────────────────────────
    draw_gene_model(axes[0], exons, view_start, view_end, strand,
                    highlight_exons=hi_exons,
                    exon_labels=locus.get("exon_labels"))
    axes[0].set_title(f"{gene_name}  ({chrom})", fontsize=9)
    axes[0].spines["bottom"].set_visible(False)
    axes[0].tick_params(bottom=False)


    # ── RNA-seq tracks (rows 1–3) ─────────────────────────────────────────────
    for i, ct in enumerate(CELL_TYPES):
        ax = axes[1 + i]
        x_kb, vals = bw_data[ct]
        draw_rna_track(ax, x_kb, vals, CELL_COLORS[ct])
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(bottom=False)

    # ── Constitutive skeleton tracks (rows 4–6) ───────────────────────────────
    for i, ct in enumerate(CELL_TYPES):
        ax = axes[4 + i]
        draw_skeleton_track(ax, skel_exons[ct], view_start, view_end, CELL_COLORS[ct])
        if i < len(CELL_TYPES) - 1:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(bottom=False)

    # ── PSI label tracks (rows 7–9) ───────────────────────────────────────────
    for i, ct in enumerate(CELL_TYPES):
        ax = axes[7 + i]
        draw_psi_track(ax, bed_dict[ct], view_start, view_end, color=CELL_COLORS[ct])
        for (hs, he), ct_psi in psi_annotations.items():
            if ct in ct_psi:
                mid_x = (hs + he) / 2 / 1e3
                ax.text(mid_x, 1.08, f"{ct_psi[ct]:.2f}",
                        ha="center", va="bottom", fontsize=6,
                        color=CELL_COLORS[ct], fontweight="bold")
        if i < len(CELL_TYPES) - 1:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(bottom=False)

    # ── Highlight spans across all 10 tracks ──────────────────────────────────
    for hs, he in hi_exons:
        for ax in axes:
            ax.axvspan(hs / 1e3, he / 1e3, color="#FFE0B2", alpha=0.30, zorder=0)

    # ── x-axis (RNA orientation: flip for minus-strand genes) ─────────────────
    if strand == "-":
        axes[0].set_xlim(view_end / 1e3, view_start / 1e3)
    else:
        axes[0].set_xlim(view_start / 1e3, view_end / 1e3)
    axes[-1].set_xlabel(f"Genomic position — {chrom} (kb)", fontsize=8)
    axes[-1].tick_params(axis="x", labelsize=7)

    # ── Section boxes, rotated labels, and legend (figure coordinates) ──────────
    fig.canvas.draw()

    # Cell-type legend anchored to upper-right of gene model panel
    p0 = axes[0].get_position()
    legend_handles = [mpatches.Patch(color=CELL_COLORS[ct], label=ct) for ct in CELL_TYPES]
    fig.legend(
        handles=legend_handles, fontsize=6.5, loc="upper right",
        bbox_to_anchor=(p0.x1, p0.y1 + 0.03),
        bbox_transform=fig.transFigure,
        framealpha=0.9, edgecolor="#CCCCCC",
        ncol=len(CELL_TYPES), handlelength=1.0, handleheight=0.8, borderpad=0.5,
    )

    sections = [
        ("RNA-seq",         1, 3),
        ("Const. skeleton", 4, 6),
        ("PSI labels",      7, 9),
    ]
    label_x_offset = 0.05  # figure-coord distance left of the axes left edge
    for sec_label, r0, r1 in sections:
        p_top = axes[r0].get_position()
        p_bot = axes[r1].get_position()
        x0 = p_top.x0
        y0 = p_bot.y0
        w  = p_top.x1 - p_top.x0
        h  = p_top.y1 - p_bot.y0

        # add extra headroom for PSI labels
        if sec_label == "PSI labels":
            extra = 0.01  # figure coords; tweak as needed
            y0 -= 0.0
            h  += extra

        rect = plt.Rectangle(
            (x0, y0), w, h,
            transform=fig.transFigure,
            fill=False, edgecolor="#BBBBBB", linewidth=0.7,
            clip_on=False,
        )
        fig.add_artist(rect)
        fig.text(
            x0 - label_x_offset, y0 + h / 2, sec_label,
            transform=fig.transFigure,
            rotation=90, va="center", ha="center",
            fontsize=7, color="#444444",
        )

    # Gene model label at same x as section labels
    fig.text(
        p0.x0 - label_x_offset, p0.y0 + (p0.y1 - p0.y0) / 2, "Gene\nmodel",
        transform=fig.transFigure,
        rotation=0, va="center", ha="center",
        fontsize=7, color="#444444",
    )

    out = os.path.join(OUT_DIR, f"panel_C_{locus_key}.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    targets = sys.argv[1:] if len(sys.argv) > 1 else list(LOCI.keys())
    for key in targets:
        if key not in LOCI:
            print(f"Unknown locus '{key}'. Available: {list(LOCI.keys())}")
            continue
        plot_locus(key)
