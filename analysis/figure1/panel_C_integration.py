"""
panel_C_integration.py

Figure 1C: rMATS PSI labels across cell types at representative alternatively-spliced loci.

Two example loci, each rendered as a 4-track genome browser (gene model + B Cell + Macrophage + T Cell):

  panel_C_ptprc.pdf  — Ptprc (CD45) cassette exons, chr1:138112000-128000
                        B Cell (high PSI) vs Macrophage (very low) shows immune-cell-specific splicing

  panel_C_pkm.pdf    — Pkm (pyruvate kinase) exon 9/10 region, chr9
                        Macrophage-specific inclusion of exon at 59675028-59675195

PSI track rendering per cell type:
  - binary skeleton from B_Cell_dominant.gtf (shared reference for all tracks)
  - SE event bins overwritten with measured mean PSI (colored by PSI colormap)
  - Constitutive-exon bins: cell-type color at PSI=1
  - Intron/intergenic bins: gray at PSI=0

Output: figure1_out/panel_C_{ptprc,pkm}.pdf
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = os.path.join(os.path.dirname(__file__), "../..")
GTF_PATH   = os.path.join(BASE, "results/stringtie_per_rep/dominant_isoforms/B_Cell_dominant.gtf")
RMATS_DIRS = {
    "B Cell":     os.path.join(BASE, "results/B_only"),
    "Macrophage": os.path.join(BASE, "results/Macrophage_only"),
    "T Cell":     os.path.join(BASE, "results/T_only"),
}
OUT_DIR = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Visual constants ───────────────────────────────────────────────────────────
CELL_TYPES = ["B Cell", "Macrophage", "T Cell"]
CELL_COLORS = {
    "B Cell":     "#2166AC",
    "Macrophage": "#4DAC26",
    "T Cell":     "#D6604D",
}
PSI_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "psi_bwr", ["#2166AC", "#F7F7F7", "#C0392B"]
)
BIN_SIZE = 32

# ── Locus definitions ─────────────────────────────────────────────────────────
LOCI = {
    "ptprc": {
        "label":      "Ptprc",
        "gene_name":  "Ptprc",
        "chrom":      "chr1",
        "strand":     "-",
        "view_start": 138_112_000,
        "view_end":   138_128_500,
        # Cassette exons to annotate; PSI varies strongly across cell types
        "highlight_exons": [
            (138_117_753, 138_117_894),
            (138_120_019, 138_120_166),
        ],
    },
    "pkm": {
        "label":      "Pkm",
        "gene_name":  "Pkm",
        "chrom":      "chr9",
        "strand":     "+",
        "view_start": 59_673_000,
        "view_end":   59_680_500,
        "highlight_exons": [
            (59_675_028, 59_675_195),   # Mac-specific exon
        ],
    },
}


# ── Data helpers ──────────────────────────────────────────────────────────────

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
    return sorted(set(exons))


def make_binary_track(exons, view_start, view_end, bin_size=BIN_SIZE):
    n_bins = (view_end - view_start) // bin_size
    track = np.zeros(n_bins, dtype=float)
    for es, ee in exons:
        b0 = max(0, (es - view_start) // bin_size)
        b1 = min(n_bins, (ee - view_start + bin_size - 1) // bin_size)
        track[b0:b1] = 1.0
    return track


def _mean_psi(s):
    if pd.isna(s) or s == "NA":
        return np.nan
    vals = [float(x) for x in str(s).split(",") if x != "NA"]
    return float(np.nanmean(vals)) if vals else np.nan


def _sum_reads(s):
    if pd.isna(s):
        return 0.0
    try:
        return sum(float(x) for x in str(s).split(",") if x not in ("", "NA"))
    except Exception:
        return 0.0


def load_se_events(rmats_dir, chrom, view_start, view_end, min_reads=20):
    path = os.path.join(rmats_dir, "SE.MATS.JC.txt")
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path, sep="\t")
    df = df[df["chr"] == chrom].copy()
    df = df[(df["exonStart_0base"] >= view_start) & (df["exonEnd"] <= view_end)]
    df["psi"] = df["IncLevel1"].apply(_mean_psi)
    df["ijc"] = df["IJC_SAMPLE_1"].apply(_sum_reads)
    df = df.dropna(subset=["psi"])
    df = df[df["ijc"] >= min_reads]
    # Keep highest-coverage event per unique exon
    df = (df.sort_values("ijc", ascending=False)
            .drop_duplicates(subset=["exonStart_0base", "exonEnd"]))
    return [(int(r["exonStart_0base"]), int(r["exonEnd"]), float(r["psi"]))
            for _, r in df.iterrows()]


def apply_psi_overwrite(binary_track, se_events, view_start, bin_size=BIN_SIZE):
    track = binary_track.copy()
    n_bins = len(track)
    for es, ee, psi in se_events:
        b0 = max(0, (es - view_start) // bin_size)
        b1 = min(n_bins, (ee - view_start + bin_size - 1) // bin_size)
        track[b0:b1] = psi
    return track


# ── Drawing helpers ───────────────────────────────────────────────────────────

def draw_gene_model(ax, exons, view_start, view_end, strand,
                    highlight_exons=(), bin_size=BIN_SIZE):
    x_lims = (view_start / 1e3, view_end / 1e3)
    ax.set_xlim(*x_lims)
    ax.set_ylim(-0.4, 1.5)
    ax.plot(x_lims, [0.5, 0.5], color="#555555", lw=0.9, zorder=1)

    n_ticks = max(1, (view_end - view_start) // 2000)
    for tp in np.linspace(view_start + 1000, view_end - 1000, n_ticks):
        dx = 400 if strand == "+" else -400
        ax.annotate("", xy=((tp + dx) / 1e3, 0.5), xytext=(tp / 1e3, 0.5),
                    arrowprops=dict(arrowstyle="-|>", color="#AAAAAA",
                                   lw=0.5, mutation_scale=6))

    hi_set = {(es, ee) for es, ee in highlight_exons}
    for es, ee in exons:
        is_hi = (es, ee) in hi_set or any(
            (hs <= es and ee <= he) for hs, he in hi_set
        )
        fc = "#D6604D" if is_hi else "#2166AC"
        ec = "#A93226" if is_hi else "#1A4E7A"
        rect = mpatches.FancyBboxPatch(
            (es / 1e3, 0.1), (ee - es) / 1e3, 0.8,
            boxstyle="square,pad=0", facecolor=fc, edgecolor=ec, lw=0.7, zorder=2,
        )
        ax.add_patch(rect)

    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)


def draw_psi_track(ax, psi_track, binary_track, view_start, view_end,
                   bin_size=BIN_SIZE, highlight_exons=()):
    n_bins = len(psi_track)
    x_kb   = view_start / 1e3 + (np.arange(n_bins) + 0.5) * bin_size / 1e3

    for xi, (xv, psi, binary) in enumerate(zip(x_kb, psi_track, binary_track)):
        if binary < 0.5 and psi < 0.01:
            c = "#E8E8E8"   # intron/intergenic
            h = 0.03        # tiny stub so intron region is visible
        elif psi < 0.01:
            c = "#CCCCCC"   # exon bin at PSI~0 (skipped)
            h = psi if psi > 0 else 0.0
        else:
            c = PSI_CMAP(psi)
            h = psi
        if h > 0:
            ax.bar(xv, h, width=bin_size / 1e3, color=c, lw=0, zorder=2)

    ax.set_xlim(view_start / 1e3, view_end / 1e3)
    ax.set_ylim(-0.05, 1.15)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "", "1"], fontsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=6)


def annotate_exons(axes_dict, se_events_dict, highlight_exons, view_start):
    """Add PSI value labels above the highlighted cassette exon(s)."""
    for hs, he in highlight_exons:
        mid_x = (hs + he) / 2 / 1e3
        for ct, ax in axes_dict.items():
            events = se_events_dict.get(ct, [])
            for es, ee, psi in events:
                if abs(es - hs) <= 50 and abs(ee - he) <= 50:
                    ax.text(mid_x, 1.08, f"{psi:.2f}",
                            ha="center", va="bottom", fontsize=6,
                            color=CELL_COLORS[ct], fontweight="bold")
                    break


# ── Main plot function ────────────────────────────────────────────────────────

def plot_locus(locus_key):
    locus = LOCI[locus_key]
    chrom      = locus["chrom"]
    gene_name  = locus["gene_name"]
    strand     = locus["strand"]
    view_start = locus["view_start"]
    view_end   = locus["view_end"]
    hi_exons   = locus.get("highlight_exons", [])

    print(f"\n── {locus['label']} ({chrom}:{view_start}-{view_end}) ──")

    exons  = parse_gtf_exons(GTF_PATH, chrom, gene_name, view_start, view_end)
    binary = make_binary_track(exons, view_start, view_end)
    print(f"  GTF exons in view: {len(exons)}")

    psi_tracks   = {}
    se_events_by = {}
    for ct in CELL_TYPES:
        evts = load_se_events(RMATS_DIRS[ct], chrom, view_start, view_end)
        se_events_by[ct] = evts
        psi_tracks[ct] = apply_psi_overwrite(binary, evts, view_start)
        psi_vals = [psi for _, _, psi in evts]
        print(f"  {ct}: {len(evts)} SE events  PSI range {min(psi_vals):.2f}-{max(psi_vals):.2f}"
              if psi_vals else f"  {ct}: 0 SE events")

    # ── Layout ────────────────────────────────────────────────────────────────
    n_tracks = 1 + len(CELL_TYPES)
    fig, axes = plt.subplots(
        n_tracks, 1, figsize=(7, 1.4 * n_tracks + 0.3),
        sharex=True,
        gridspec_kw={"height_ratios": [0.9] + [1.0] * len(CELL_TYPES)},
    )
    fig.subplots_adjust(hspace=0.04)

    # ── Gene model ────────────────────────────────────────────────────────────
    draw_gene_model(axes[0], exons, view_start, view_end, strand,
                    highlight_exons=hi_exons)
    axes[0].set_title(
        f"{gene_name} — PSI labels by cell type  ({chrom})", fontsize=9
    )
    axes[0].set_ylabel("Gene\nmodel", fontsize=7, rotation=0, labelpad=46, va="center")

    # ── PSI tracks ────────────────────────────────────────────────────────────
    ct_axes = {}
    for ti, ct in enumerate(CELL_TYPES):
        ax = axes[ti + 1]
        ct_axes[ct] = ax
        draw_psi_track(ax, psi_tracks[ct], binary,
                       view_start, view_end, highlight_exons=hi_exons)
        ax.set_ylabel(ct, fontsize=7, rotation=0, labelpad=46, va="center",
                      color=CELL_COLORS[ct], fontweight="bold")
        if ti < len(CELL_TYPES) - 1:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(bottom=False)

    annotate_exons(ct_axes, se_events_by, hi_exons, view_start)

    # Highlight cassette exon spans with vertical shading on all tracks
    for hs, he in hi_exons:
        for ax in axes:
            ax.axvspan(hs / 1e3, he / 1e3, color="#FFE0B2", alpha=0.30, zorder=0)

    # ── x-axis ────────────────────────────────────────────────────────────────
    axes[-1].set_xlabel(f"Genomic position — {chrom} (kb)", fontsize=8)
    axes[-1].tick_params(axis="x", labelsize=7)

    # ── PSI colorbar ──────────────────────────────────────────────────────────
    sm = plt.cm.ScalarMappable(cmap=PSI_CMAP, norm=mcolors.Normalize(0, 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[1:], orientation="vertical",
                        fraction=0.012, pad=0.01, shrink=0.8)
    cbar.set_label("PSI", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

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
