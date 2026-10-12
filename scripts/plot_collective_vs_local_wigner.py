#!/usr/bin/env python3
"""Figure: local vs collective edge-mode Wigner functions of the final state.

Left panel:  single-site (site 0) Wigner function of the trajectory-averaged
             final state (topological sector, edge-selective loss, d=0, sample 0)
             -- the production observable, with a genuine negative region.
Right panel: collective edge-mode Wigner function of the same state after the
             validated mode rotation -- positive everywhere (no negativity).

This visualizes the null result of the collective-mode analysis: the finite-time
non-Gaussianity lives in the site basis and does not survive rotation into the
edge-mode basis.

The plotted states are read from data/collective_wigner_figure_state.json
(written by run_collective_wigner.py --dump-figure-state); this script runs
no dynamics and imports no producer module.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.style.use(["seaborn-v0_8-ticks", "seaborn-v0_8-colorblind"])
import numpy as np
import qutip as qt
from matplotlib import patheffects as path_effects

ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = ROOT / "data" / "collective_wigner_figure_state.json"
OUT_PDF = ROOT / "figures" / "collective_vs_local_wigner.pdf"
OUT_PNG = ROOT / "figures" / "collective_vs_local_wigner.png"

EXPECTED_CONFIG = {
    "n_cells": 2, "n_max": 2, "sector": "chiral_topological",
    "loss_profile": "edge_low_loss", "disorder_level": 0.0,
    "sample": 0, "draw_seed": 20260814,
    "trajectory_seed": 20260914, "ntraj": 10,
}


def load_block(payload: dict, key: str) -> qt.Qobj:
    """Rebuild a 3x3 density matrix from its archived real/imag parts.

    Args:
        payload (dict): The figure-state artifact.
        key (str): Block entry name.

    Returns:
        qutip.Qobj: The reduced density matrix (dimensionless entries).
    """
    block = payload[key]
    mat = np.array(block["real"]) + 1j * np.array(block["imag"])
    return qt.Qobj(mat)


def main():
    payload = json.loads(STATE_PATH.read_text())
    config = payload["config"]
    for field, want in EXPECTED_CONFIG.items():
        got = config[field]
        assert got == want, f"figure-state config drift: {field} {got} != {want}"
    rho_local = load_block(payload, "rho_local_site0")
    rho_edge = load_block(payload, "rho_collective_edge")

    grid = np.linspace(-4.0, 4.0, 161)
    w_local = qt.wigner(rho_local, grid, grid)
    w_edge = qt.wigner(rho_edge, grid, grid)

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.9), constrained_layout=True)
    for ax, w, letter in [
        (axes[0], w_local, "(a)"),
        (axes[1], w_edge, "(b)"),
    ]:
        lev = np.linspace(w.min(), w.max(), 41)
        cf = ax.contourf(grid, grid, w, levels=lev, cmap="RdBu_r")
        ax.contour(grid, grid, w, levels=[0.0], colors="k", linewidths=0.7)
        ax.set_xlabel(r"$\mathrm{Re}\,\alpha$")
        ax.set_ylabel(r"$\mathrm{Im}\,\alpha$")
        ax.text(0.02, 0.98, letter, transform=ax.transAxes,
                va="top", ha="left", fontsize=11, color="white",
                path_effects=[path_effects.withStroke(linewidth=1.5, foreground="black")])
        fig.colorbar(cf, ax=ax)
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    print(f"saved {OUT_PDF}")
    print(f"saved {OUT_PNG}")
    print(f"local min W = {w_local.min():.4f}, edge min W = {w_edge.min():.2e}")


if __name__ == "__main__":
    main()
