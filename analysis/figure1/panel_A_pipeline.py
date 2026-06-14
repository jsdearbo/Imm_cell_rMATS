"""
panel_A_pipeline.py

Figure 1A: Cell-type-specific splicing label generation pipeline flowchart.

Two output files:
  figure1_out/panel_A_pipeline_multitask.pdf   -- stops at chromosome-split NPZ shards
  figure1_out/panel_A_pipeline_singletask.pdf  -- adds single-task NPZ extraction step
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

OUT_DIR = os.path.join(os.path.dirname(__file__), "figure1_out")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Palette ───────────────────────────────────────────────────────────────────
C_IN_FC = "#78909C"; C_IN_EC = "#37474F"   # input (blue-gray)
C_BL_FC = "#1565C0"; C_BL_EC = "#0D47A1"  # annotation branch (blue)
C_OR_FC = "#E64A19"; C_OR_EC = "#BF360C"  # rMATS branch (orange)
C_PU_FC = "#6A1B9A"; C_PU_EC = "#4A148C"  # soft-label synthesis (purple, emphasized)
C_GR_FC = "#2E7D32"; C_GR_EC = "#1B5E20"  # NPZ output (green)
C_G2_FC = "#558B2F"; C_G2_EC = "#33691E"  # single-task output (olive-green)
C_ARR   = "#424242"

CT_COLORS = {
    "B Cell":     "#2166AC",
    "T Cell":     "#D6604D",
    "Macrophage": "#4DAC26",
}

MAIN_W = 8.8   # full-width boxes
BRAN_W = 4.1   # branch boxes


# ── Drawing primitives ────────────────────────────────────────────────────────

def fbox(ax, cx, cy, w, h, txt, fc, ec,
         fs=7.5, bold=False, tc="white", lw=1.1, zorder=2):
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.03",
        facecolor=fc, edgecolor=ec, linewidth=lw,
        clip_on=False, zorder=zorder,
    ))
    ax.text(cx, cy, txt, ha="center", va="center",
            fontsize=fs, color=tc,
            fontweight="bold" if bold else "normal",
            linespacing=1.45, zorder=zorder + 1)


def vline(ax, x, y0, y1, color=C_ARR, lw=1.0):
    ax.plot([x, x], [y0, y1], color=color, lw=lw,
            solid_capstyle="butt", zorder=1)


def hline(ax, x0, x1, y, color=C_ARR, lw=1.0):
    ax.plot([x0, x1], [y, y], color=color, lw=lw,
            solid_capstyle="butt", zorder=1)


def arrow_down(ax, x, y_from, y_to, color=C_ARR):
    ax.annotate("", xy=(x, y_to), xytext=(x, y_from),
                arrowprops=dict(arrowstyle="-|>", color=color,
                                lw=1.0, mutation_scale=9),
                zorder=5)


# ── Tensor schematic ──────────────────────────────────────────────────────────

def tensor_schematic(ax, y_top, width=MAIN_W, height=2.5):
    """Draw the NPZ tensor schematic below the flowchart boxes."""
    cx = 5.0
    x0 = cx - width / 2
    x1 = cx + width / 2
    y0 = y_top - height

    # Outer border
    ax.add_patch(FancyBboxPatch(
        (x0, y0), width, height,
        boxstyle="round,pad=0.03",
        facecolor="#F5F5F5", edgecolor="#AAAAAA", linewidth=0.8,
        clip_on=False, zorder=1,
    ))
    ax.text(cx, y_top - 0.15, "Model-facing tensor (NPZ format)",
            ha="center", va="top", fontsize=7, color="#333333",
            fontweight="bold", zorder=3)

    # Track layout
    LABEL_W = 1.75
    track_x0 = x0 + LABEL_W + 0.08
    track_w = width - LABEL_W - 0.22
    N = 32
    bw = track_w / N

    n_rows = 5
    row_gap = 0.10
    row_h_avail = height - 0.38 - row_gap * (n_rows - 1)
    rh = row_h_avail / n_rows

    # Demo data: 2 constitutive exons flanking 1 alt exon
    # Bins: 0-1 intergenic, 2-7 const exon1, 8-11 intron, 12-17 alt exon,
    #        18-21 intron, 22-27 const exon2, 28-31 intergenic
    np.random.seed(42)
    const1 = np.zeros(N)
    const1[2:8] = 1.0; const1[22:28] = 1.0

    def _make_psi(alt_val, noise_scale=0.04):
        arr = const1.copy()
        arr[12:18] = alt_val
        arr += np.random.rand(N) * noise_scale
        return np.clip(arr, 0.0, 1.0)

    seq  = np.clip(0.45 + np.random.rand(N) * 0.35, 0, 1)
    alt_B = _make_psi(0.10)
    alt_T = _make_psi(0.78)
    alt_M = _make_psi(0.43)
    mask  = np.zeros(N)
    mask[2:8] = 1; mask[12:18] = 1; mask[22:28] = 1

    rows = [
        ("Sequence (one-hot)", "#37474F", seq),
        ("B cell PSI",         CT_COLORS["B Cell"],     alt_B),
        ("T cell PSI",         CT_COLORS["T Cell"],     alt_T),
        ("Macrophage PSI",     CT_COLORS["Macrophage"], alt_M),
        ("Mask",               "#78909C",               mask),
    ]

    for ri, (label, color, vals) in enumerate(rows):
        ry = y_top - 0.32 - ri * (rh + row_gap) - rh / 2
        ax.text(x0 + LABEL_W - 0.06, ry, label,
                ha="right", va="center", fontsize=5.8, color="#333333",
                zorder=3)
        for bi, v in enumerate(vals):
            bx = track_x0 + bi * bw
            is_signal = (label == "Sequence (one-hot)") or (v > 0.07)
            clr = color if is_signal else "#E0E0E0"
            alpha = float(np.clip(v + 0.18, 0.18, 1.0))
            ax.add_patch(mpatches.Rectangle(
                (bx, ry - rh * 0.42), bw * 0.84, rh * 0.84,
                facecolor=clr, edgecolor="none", alpha=alpha,
                clip_on=True, zorder=2,
            ))

    ax.text(x1 - 0.1, y0 + 0.08,
            "32-bp bins  ·  524 kb window  ·  3 cell-type tasks",
            ha="right", va="bottom", fontsize=5.5,
            color="#666666", fontstyle="italic", zorder=3)


# ── Main figure function ──────────────────────────────────────────────────────

def make_figure(singletask: bool = False):
    suffix = "singletask" if singletask else "multitask"

    fig_h = 10.5 if singletask else 9.5
    y_lo  = 2.8  if singletask else 3.8

    fig, ax = plt.subplots(figsize=(6.5, fig_h))
    ax.set_xlim(0, 10)
    ax.set_ylim(y_lo, 17.3)
    ax.axis("off")

    # ── Input box ──────────────────────────────────────────────────────────────
    in_cy, in_h = 16.2, 1.05
    fbox(ax, 5, in_cy, MAIN_W, in_h,
         "RNA-seq BAMs\nB cells  ·  T cells  ·  Macrophages\nbiological replicates per cell type",
         C_IN_FC, C_IN_EC, fs=8.0, bold=True)
    in_bot = in_cy - in_h / 2  # 15.675

    # ── Fork ───────────────────────────────────────────────────────────────────
    FORK_Y = in_bot - 0.42      # 15.255 — horizontal bar
    vline(ax, 5, in_bot, FORK_Y)
    hline(ax, 2.5, 7.5, FORK_Y)

    # Branch column headers (colored, between fork bar and first boxes)
    ax.text(2.5, FORK_Y - 0.12, "Transcript annotation",
            ha="center", va="top", fontsize=7.0,
            color=C_BL_EC, fontweight="bold", zorder=6)
    ax.text(7.5, FORK_Y - 0.12, "Splicing quantification",
            ha="center", va="top", fontsize=7.0,
            color=C_OR_EC, fontweight="bold", zorder=6)

    # ── Branch box 1 ───────────────────────────────────────────────────────────
    B1_H = 0.88
    B1_TOP = FORK_Y - 0.72       # 14.535
    B1_CY  = B1_TOP - B1_H / 2  # 14.095

    arrow_down(ax, 2.5, FORK_Y, B1_TOP)
    arrow_down(ax, 7.5, FORK_Y, B1_TOP)

    fbox(ax, 2.5, B1_CY, BRAN_W, B1_H,
         "StringTie transcript assembly\n(per replicate, per cell type)",
         C_BL_FC, C_BL_EC)
    fbox(ax, 7.5, B1_CY, BRAN_W, B1_H,
         "rMATS splice event detection\n(per cell type)",
         C_OR_FC, C_OR_EC)
    B1_BOT = B1_CY - B1_H / 2   # 13.655

    # ── Branch box 2 ───────────────────────────────────────────────────────────
    B2_H   = 1.12
    B2_GAP = 0.45
    B2_TOP = B1_BOT - B2_GAP     # 13.205
    B2_CY  = B2_TOP - B2_H / 2  # 12.645

    arrow_down(ax, 2.5, B1_BOT, B2_TOP)
    arrow_down(ax, 7.5, B1_BOT, B2_TOP)

    fbox(ax, 2.5, B2_CY, BRAN_W, B2_H,
         "Dominant isoform selection\n→ constitutive skeleton\nexon = 1  ·  intron = 0",
         C_BL_FC, C_BL_EC)
    fbox(ax, 7.5, B2_CY, BRAN_W, B2_H,
         "SE · A5SS · A3SS · RI · MXE\nper-sample PSI intervals\n+ read-support confidence\n(low-support events masked)",
         C_OR_FC, C_OR_EC)
    B2_BOT = B2_CY - B2_H / 2   # 12.085

    # ── Join ───────────────────────────────────────────────────────────────────
    JOIN_Y = B2_BOT - 0.35       # 11.735
    vline(ax, 2.5, B2_BOT, JOIN_Y)
    vline(ax, 7.5, B2_BOT, JOIN_Y)
    hline(ax, 2.5, 7.5, JOIN_Y)

    # ── Soft-label synthesis (emphasized) ──────────────────────────────────────
    SYN_H   = 1.50
    SYN_GAP = 0.45
    SYN_TOP = JOIN_Y - SYN_GAP   # 11.285
    SYN_CY  = SYN_TOP - SYN_H / 2  # 10.535

    arrow_down(ax, 5, JOIN_Y, SYN_TOP)
    fbox(ax, 5, SYN_CY, MAIN_W, SYN_H,
         "Cell-type-specific soft-label synthesis\n"
         "PSI values overwrite corresponding skeleton regions\n"
         "→ replicate aggregation  →  master BED per cell type",
         C_PU_FC, C_PU_EC, fs=8.0, bold=True, lw=1.6)
    SYN_BOT = SYN_CY - SYN_H / 2  # 9.785

    # ── NPZ conversion ─────────────────────────────────────────────────────────
    NPZ_H   = 1.12
    NPZ_GAP = 0.45
    NPZ_TOP = SYN_BOT - NPZ_GAP   # 9.335
    NPZ_CY  = NPZ_TOP - NPZ_H / 2  # 8.775

    arrow_down(ax, 5, SYN_BOT, NPZ_TOP)
    fbox(ax, 5, NPZ_CY, MAIN_W, NPZ_H,
         "NPZ tensor conversion\n"
         "sequence + PSI + mask  ·  32-bp bins  ·  3 cell-type tasks\n"
         "chromosome-held-out  train / validation / test",
         C_GR_FC, C_GR_EC, fs=8.0)
    NPZ_BOT = NPZ_CY - NPZ_H / 2  # 8.215

    # ── Single-task extraction (optional) ──────────────────────────────────────
    if singletask:
        SNG_H   = 0.88
        SNG_GAP = 0.45
        SNG_TOP = NPZ_BOT - SNG_GAP   # 7.765
        SNG_CY  = SNG_TOP - SNG_H / 2  # 7.325

        arrow_down(ax, 5, NPZ_BOT, SNG_TOP)
        fbox(ax, 5, SNG_CY, MAIN_W, SNG_H,
             "Single-task NPZ extraction\n"
             "B cell-only  ·  T cell-only  ·  macrophage-only labels",
             C_G2_FC, C_G2_EC, fs=8.0)
        SNG_BOT = SNG_CY - SNG_H / 2  # 6.885
        tensor_top = SNG_BOT - 0.50   # 6.385
    else:
        tensor_top = NPZ_BOT - 0.50   # 7.715

    # ── Tensor schematic ───────────────────────────────────────────────────────
    tensor_schematic(ax, tensor_top, width=MAIN_W, height=2.5)

    # ── Save ───────────────────────────────────────────────────────────────────
    out = os.path.join(OUT_DIR, f"panel_A_pipeline_{suffix}.pdf")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    make_figure(singletask=False)
    make_figure(singletask=True)
