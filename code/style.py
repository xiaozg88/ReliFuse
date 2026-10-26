# -*- coding: utf-8 -*-
"""Shared plotting style for v46 figures (SCI Q1 look, Okabe-Ito palette)."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.normpath(os.path.join(HERE, "..", "figures"))
DATADIR = os.path.normpath(os.path.join(HERE, "..", "data"))
os.makedirs(FIGDIR, exist_ok=True)

# Okabe-Ito colour-blind safe palette
C = {
    "blue": "#0072B2",
    "sky": "#56B4E9",
    "green": "#009E73",
    "amber": "#E69F00",
    "orange": "#D55E00",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "grey": "#7F7F7F",
    "lightgrey": "#ECECEC",
    "pale": "#F5F5F5",
    "ink": "#222222",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "axes.titleweight": "bold",
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "figure.dpi": 120,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "axes.grid": True,
    "grid.alpha": 0.18,
    "grid.linewidth": 0.5,
    "axes.axisbelow": True,
    "lines.linewidth": 1.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save(fig, name):
    """Save a figure as both PDF (vector) and PNG (400 dpi)."""
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"),
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved", name)


def rbox(ax, x, y, w, h, text="", fc="white", ec=C["ink"], lw=0.9,
         dashed=False, fs=8, bold=False, tc=None, radius=0.018, zorder=3,
         ha="center"):
    """Rounded rectangle with centred (or left) text; coordinates in axes units."""
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0,rounding_size={radius}",
        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=zorder,
        linestyle=(0, (3, 2)) if dashed else "solid",
    )
    ax.add_patch(box)
    if text:
        tx = x + w / 2 if ha == "center" else x + 0.012
        ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fs,
                fontweight="bold" if bold else "normal",
                color=tc or C["ink"], zorder=zorder + 1, linespacing=1.35)
    return box


def arrow(ax, p0, p1, color=None, lw=1.0, ls="solid", zorder=2, rad=0.0):
    a = FancyArrowPatch(
        p0, p1, arrowstyle="-|>", mutation_scale=8, shrinkA=1.5, shrinkB=1.5,
        linewidth=lw, color=color or C["grey"], linestyle=ls, zorder=zorder,
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(a)
    return a


def blank_axes(fig, xlim=(0, 1), ylim=(0, 1)):
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis("off")
    return ax


def panel_label(ax, s, x=-0.08, y=1.06, fs=10):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=fs,
            fontweight="bold", va="top", ha="left", color=C["ink"])
