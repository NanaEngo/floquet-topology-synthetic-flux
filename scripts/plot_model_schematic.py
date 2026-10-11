#!/usr/bin/env python3
"""Model schematic (roadmap G11, Figure 1).

Three panels: (a) the driven chiral SSH chain with edge/bulk loss channels and
the local observable, (b) the two-step Floquet drive protocol, (c) the even-cat
edge preparation whose local Wigner negativity is the tracked resource.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
plt.style.use(["seaborn-v0_8-ticks", "seaborn-v0_8-colorblind"])
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT_PDF = ROOT / "figures" / "model_schematic.pdf"
OUT_PNG = ROOT / "figures" / "model_schematic.png"


def draw_chain(ax):
    """Two-cell (four-mode) chain: solid = intracell J_A, dashed = intercell J_B."""
    n_cells = 2
    sites = []
    for c in range(n_cells):
        sites.append((c * 2.0, 0.0))          # A site
        sites.append((c * 2.0 + 1.0, 0.0))    # B site
    for i, (x, y) in enumerate(sites):
        ax.plot(x, y, "o", ms=16, color="#1f77b4", zorder=3)
        ax.text(x, y - 0.33, f"{i}", ha="center", va="top", fontsize=9)
    # bonds: intracell (A-B within cell) solid; intercell (A of next) dashed
    for c in range(n_cells):
        ax.plot([c * 2.0, c * 2.0 + 1.0], [0, 0], "-", lw=2.2, color="#d62728",
                zorder=1)
    for c in range(n_cells - 1):
        ax.plot([c * 2.0 + 1.0, (c + 1) * 2.0], [0, 0], "--", lw=2.2, color="#0072B2",
                zorder=1)
    ax.text(0.5, 0.30, r"$J_A$ (step A)", ha="center", fontsize=10, color="#d62728")
    ax.text(1.5, 0.30, r"$J_B$ (step B)", ha="center", fontsize=10, color="#0072B2")
    # loss arrows: weak at edges, strong in bulk
    for x, tag, gamma in [(0.0, "edge", r"$\gamma_{\rm edge}$"),
                          (1.0, "bulk", r"$\gamma_{\rm bulk}$"),
                          (2.0, "bulk", r"$\gamma_{\rm bulk}$"),
                          (3.0, "edge", r"$\gamma_{\rm edge}$")]:
        ax.annotate("", xy=(x, -0.95), xytext=(x, -0.28),
                    arrowprops=dict(arrowstyle="-|>", lw=1.4,
                                    color="#7f7f7f" if tag == "bulk" else "#ff7f0e"))
        ax.text(x, -1.15, gamma, ha="center", fontsize=9,
                color="#7f7f7f" if tag == "bulk" else "#ff7f0e")
    # observable marker at site 0
    ax.annotate("", xy=(0.0, 0.85), xytext=(0.0, 0.30),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#9467bd"))
    ax.text(0.0, 1.05, r"$N_-(\rho_0)$ local", ha="center", fontsize=10, color="#9467bd")
    ax.set_xlim(-0.9, 3.9)
    ax.set_ylim(-1.5, 1.6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    ax.text(0.02, 0.98, "(a)", transform=ax.transAxes,
            va="top", ha="left", fontsize=12)


def draw_drive(ax):
    """Two-step drive: step A (intracell) then step B (intercell), period T = 2*pi."""
    T = 2.0 * np.pi
    t = np.linspace(0.0, 2 * T, 1000)
    step = np.where((t % T) < T / 2, 1.0, 0.0)      # step A active
    ax.plot(t / T, step, lw=2.0, color="#d62728", label="step A (intracell $J_A$)")
    ax.plot(t / T, 1.0 - step, lw=2.0, color="#0072B2", label="step B (intercell $J_B$)")
    ax.axvline(0.5, color="k", lw=0.8, ls=":")
    ax.axvline(1.0, color="k", lw=0.8, ls=":")
    ax.axvline(1.5, color="k", lw=0.8, ls=":")
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.08, 1.08)
    ax.set_xlabel(r"$t/T$")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["off", "on"])
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    ax.text(0.02, 0.98, "(b)", transform=ax.transAxes,
            va="top", ha="left", fontsize=12)


def draw_cat(ax):
    """Even-cat edge preparation in phase space (schematic)."""
    alpha = 1.0
    x = np.linspace(-3.2, 3.2, 400)
    y = np.linspace(-3.2, 3.2, 400)
    X, Y = np.meshgrid(x, y)
    lobe = np.exp(-((X - alpha) ** 2 + Y ** 2)) + np.exp(-((X + alpha) ** 2 + Y ** 2))
    fringe = 2 * np.exp(-(X ** 2 + Y ** 2)) * np.cos(2 * alpha * Y)
    w = lobe + fringe
    ax.contourf(X, Y, w, levels=30, cmap="RdBu")
    ax.contour(X, Y, w, levels=[0.0], colors="k", lw=1.2)
    ax.set_xlabel(r"$\mathrm{Re}\,\alpha$")
    ax.set_ylabel(r"$\mathrm{Im}\,\alpha$")
    ax.set_xticks([-2, 0, 2])
    ax.set_yticks([-2, 0, 2])
    ax.set_aspect("equal")
    ax.text(0.02, 0.98, "(c)", transform=ax.transAxes,
            va="top", ha="left", fontsize=12)


def main():
    fig = plt.figure(figsize=(12.0, 3.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0])
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2])
    draw_chain(ax_a)
    draw_drive(ax_b)
    draw_cat(ax_c)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"Saved {OUT_PDF}")
    print(f"Saved {OUT_PNG}")


if __name__ == "__main__":
    main()
