# -*- coding: utf-8 -*-
"""Figure 1 (v46): ReliFuse architecture schematic, SCI Q1 style."""

import matplotlib.pyplot as plt

from style import C, rbox, arrow, save, blank_axes


def build():
    fig = plt.figure(figsize=(7.4, 5.9))
    ax = blank_axes(fig, xlim=(0, 1), ylim=(0, 1))

    # ------------------------------------------------------------------ title
    ax.text(0.5, 0.985,
            "ReliFuse: multi-kernel message passing with calibrated logit fusion",
            ha="center", va="top", fontsize=9.5, fontweight="bold", color=C["ink"])

    # ------------------------------------------------------- deployment note
    rbox(ax, 0.03, 0.005, 0.94, 0.052,
         "Deployment rule — structure logits are used only when AlphaFold and "
         "ESMFold agree in sign; otherwise the pipeline falls back to the sequence branch",
         fc="white", ec=C["grey"], dashed=True, fs=6.6, tc="#444444")

    # ----------------------------------------------------------------- inputs
    rbox(ax, 0.03, 0.077, 0.44, 0.072, "Protein sequence",
         fc="#EAF6FB", ec=C["sky"], fs=8.5, bold=True)
    rbox(ax, 0.53, 0.077, 0.44, 0.072, "AlphaFold structure",
         fc="#EAF6FB", ec=C["sky"], fs=8.5, bold=True)

    # -------------------------------------------------------- sequence branch
    rbox(ax, 0.03, 0.185, 0.44, 0.355, "", fc="#F4FBF8", ec=C["green"], lw=1.0)
    rbox(ax, 0.03, 0.488, 0.44, 0.052, "Sequence branch", fc=C["green"],
         ec=C["green"], fs=8.5, bold=True, tc="white")
    rbox(ax, 0.05, 0.405, 0.40, 0.062,
         "Frozen protein language model\nESM-2 150M / 650M · ESM-C 300M",
         fc="white", ec=C["green"], dashed=True, fs=7.2)
    rbox(ax, 0.05, 0.325, 0.40, 0.055,
         r"mean-pool  $h_p=\mathrm{mean}_t(e_t)$",
         fc="white", ec=C["green"], fs=7.2)
    rbox(ax, 0.05, 0.245, 0.40, 0.055,
         r"trainable MLP head  $\rightarrow\ z_2$",
         fc="#DFF3EB", ec=C["green"], fs=7.2)
    arrow(ax, (0.25, 0.405), (0.25, 0.382), C["green"], lw=0.9)
    arrow(ax, (0.25, 0.325), (0.25, 0.302), C["green"], lw=0.9)
    ax.text(0.25, 0.208, r"$z_2$  (sequence logit)", ha="center", va="center",
            fontsize=7.6, color=C["green"], fontweight="bold")

    # ------------------------------------------------------- structure branch
    rbox(ax, 0.53, 0.185, 0.44, 0.355, "", fc="#F1F7FC", ec=C["blue"], lw=1.0)
    rbox(ax, 0.53, 0.488, 0.44, 0.052, "Structure branch", fc=C["blue"],
         ec=C["blue"], fs=8.5, bold=True, tc="white")
    rbox(ax, 0.55, 0.434, 0.40, 0.046,
         r"Residue contact graph  $G=(V,E,x)$   (C$\alpha$–C$\alpha$ < 8.0 Å,  $x_i\in\mathbb{R}^{32}$)",
         fc="white", ec=C["blue"], fs=6.8)
    ax.text(0.55, 0.418, "4 mechanism kernels per edge",
            fontsize=6.6, ha="left", va="center", color=C["ink"])
    kernels = [r"$\kappa_1$ cosine", r"$\kappa_2$ chain", r"$\kappa_3$ degree",
               r"$\kappa_4$ RBF"]
    for i, k in enumerate(kernels):
        rbox(ax, 0.55 + i * 0.101, 0.372, 0.094, 0.038, k,
             fc="#E3F0F9", ec=C["blue"], fs=6.0)
    rbox(ax, 0.55, 0.294, 0.40, 0.058,
         r"per-layer protein-level gate $\alpha^{(l)}=\mathrm{softmax}(\mathrm{MLP}_l(\bar{x})/\tau_l)$" "\n"
         r"edge weight  $g_e^{(l)}=\beta\sum_k \alpha_{l,k}\kappa_k(e)+(1-\beta)$",
         fc="#E9F3FA", ec=C["blue"], fs=6.2)
    arrow(ax, (0.63, 0.434), (0.63, 0.412), C["blue"], lw=0.9)
    arrow(ax, (0.63, 0.372), (0.63, 0.354), C["blue"], lw=0.9)
    arrow(ax, (0.63, 0.294), (0.63, 0.278), C["blue"], lw=0.9)
    rbox(ax, 0.55, 0.232, 0.40, 0.044,
         "3 × [ GraphConv + graph SE ]",
         fc="#E3F0F9", ec=C["blue"], fs=6.6)
    ax.text(0.75, 0.215, r"readout MLP  $\rightarrow\ z_1$   (structure logit)",
            ha="center", va="center", fontsize=6.8, color=C["blue"], fontweight="bold")

    # ---------------------------------------------------------- fusion module
    rbox(ax, 0.06, 0.565, 0.88, 0.272, "", fc="#FFF8EC", ec=C["amber"], lw=1.1)
    rbox(ax, 0.06, 0.790, 0.88, 0.047, "Calibrated logit fusion",
         fc=C["amber"], ec=C["amber"], fs=8.5, bold=True, tc="white")
    rbox(ax, 0.085, 0.672, 0.16, 0.098,
         "Temperature\nscaling\n" r"$z_{1c},\, z_{2c}$",
         fc="white", ec=C["amber"], fs=6.8)
    rbox(ax, 0.265, 0.672, 0.245, 0.098,
         "Gate features  (16-d)\n" 
         r"$z_{1c}, z_{2c}, p_1, p_2$, entropy, agreement," "\n"
         "disagreement, margins, pLDDT, log-length, density",
         fc="white", ec=C["amber"], fs=6.2)
    rbox(ax, 0.53, 0.672, 0.16, 0.098, "Shared trunk\nMLP 16 → 64 → 32",
         fc="white", ec=C["amber"], fs=6.8)
    rbox(ax, 0.71, 0.722, 0.20, 0.048, r"$g=\sigma(\mathrm{MLP}_g(h))$",
         fc="#FDEBD2", ec=C["amber"], fs=6.8)
    rbox(ax, 0.71, 0.668, 0.20, 0.048, r"$b=\mathrm{MLP}_b(h)$",
         fc="#FDEBD2", ec=C["amber"], fs=6.8)
    arrow(ax, (0.245, 0.721), (0.263, 0.721), C["amber"], lw=0.9)
    arrow(ax, (0.510, 0.721), (0.528, 0.721), C["amber"], lw=0.9)
    arrow(ax, (0.690, 0.746), (0.708, 0.746), C["amber"], lw=0.9)
    arrow(ax, (0.690, 0.696), (0.708, 0.692), C["amber"], lw=0.9)
    rbox(ax, 0.085, 0.582, 0.72, 0.070,
         r"$z_{\mathrm{fuse}}=z_{1c}+g\,(z_{2c}-z_{1c})+b$",
         fc="#FFFDF7", ec=C["amber"], fs=9.5, bold=True)
    rbox(ax, 0.820, 0.582, 0.095, 0.070, "≈3.2 K\ntrainable\nparameters",
         fc="#FDF3E3", ec=C["orange"], fs=6.2, tc=C["orange"], bold=True)

    # ------------------------------------------------------------------ output
    tiles = ["σ($z_{fuse}$) → p(BALF-detected)",
             "Top-K candidate ranking\n(PRM/SRM shortlist)",
             "Brier / ECE\ncalibration",
             "Decision curve\nanalysis"]
    xs = [0.06, 0.285, 0.51, 0.735]
    for x, t in zip(xs, tiles):
        rbox(ax, x, 0.875, 0.205, 0.082, t, fc="#F2F2F2", ec=C["grey"], fs=6.8)

    # --------------------------------------------------------- band connectors
    arrow(ax, (0.25, 0.149), (0.25, 0.183), C["grey"], lw=1.0)
    arrow(ax, (0.75, 0.149), (0.75, 0.183), C["grey"], lw=1.0)
    arrow(ax, (0.25, 0.540), (0.25, 0.563), C["grey"], lw=1.2)
    arrow(ax, (0.75, 0.540), (0.75, 0.563), C["grey"], lw=1.2)
    arrow(ax, (0.50, 0.837), (0.50, 0.873), C["grey"], lw=1.4)

    save(fig, "fig1_architecture")


if __name__ == "__main__":
    build()
