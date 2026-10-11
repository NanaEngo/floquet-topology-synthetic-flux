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
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.style.use(["seaborn-v0_8-ticks", "seaborn-v0_8-colorblind"])
import numpy as np
import qutip as qt
from matplotlib import patheffects as path_effects

from run_chiral_floquet_dynamics import (
    CANDIDATES,
    PERIOD,
    build_ops,
    driven_hamiltonian,
    edge_mode_vector,
    initial_state,
    rate_vector,
)
from run_collective_wigner import complete_mode_rotation, many_body_rotation, rotate_state

ROOT = Path(__file__).resolve().parents[1]
OUT_PDF = ROOT / "figures" / "collective_vs_local_wigner.pdf"
OUT_PNG = ROOT / "figures" / "collective_vs_local_wigner.png"


def main():
    n_cells, n_max, n_modes = 2, 2, 4
    sector, profile, disorder = "chiral_topological", "edge_low_loss", 0.0
    draw_seed = 20260814  # sample 0, matched with the pilot
    ops = build_ops(n_modes, n_max)
    psi0 = initial_state(n_cells, n_max, "edge_site")
    rates = rate_vector(n_modes, profile, disorder, draw_seed, 1e-4, 0.2)
    c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
    result = qt.mcsolve(
        driven_hamiltonian(ops, CANDIDATES[sector], 0.5), psi0,
        np.linspace(0.0, PERIOD, 9), c_ops, ntraj=10, seeds=draw_seed + 100,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    state = result.states[-1]

    rho_local = state.ptrace(0)
    c = edge_mode_vector(n_cells, CANDIDATES[sector])
    u = complete_mode_rotation(c)
    T = many_body_rotation(u, n_modes, n_max)
    rho_rot = rotate_state(state, T, n_modes, n_max)
    rho_edge = rho_rot.ptrace(0)

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
