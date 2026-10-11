#!/usr/bin/env python3
"""Resolution and truncation checks for the collective-mode Wigner result.

The dynamics pilot reports no resolved collective edge-mode Wigner negativity above the numerical floor
in all matched conditions, including where the local site-0 negativity is
positive.  Before that statement enters any document, rule out two artefacts:

  1. grid resolution: the 31x31 grid on [-4,4]^2 may miss a small negative
     volume; re-evaluate on a 61x61 grid and report the Wigner minimum.
  2. Fock truncation: the zero could be an n_max=2 artefact; re-evaluate at
     n_max=3 with the same representative condition.

The representative condition is the topological sector, edge_low_loss profile,
zero disorder, sample 0 (draw_seed 20260814) -- the first production record.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import qutip as qt

from run_chiral_floquet_dynamics import (
    CANDIDATES,
    PERIOD,
    build_ops,
    driven_hamiltonian,
    edge_mode_vector,
    initial_state,
    rate_vector,
)
from run_collective_wigner import (
    complete_mode_rotation,
    many_body_rotation,
    rotate_state,
    wigner_negativity,
)

ROOT = Path(__file__).resolve().parents[1]
ALPHA = 1.0


def wigner_stats(rho, grid_points):
    grid = np.linspace(-4.0, 4.0, grid_points)
    w = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    neg = float(np.sum(np.abs(w[w < 0])) * dx * dx)
    return neg, float(w.min()), float(w.max())


def run_condition(n_cells, n_max, sector, profile, disorder, ntraj, seed, t_stop=PERIOD):
    n_modes = 2 * n_cells
    ops = build_ops(n_modes, n_max)
    psi0 = initial_state(n_cells, n_max, "edge_site")
    rates = rate_vector(n_modes, profile, disorder, seed, 1e-4, 0.2)
    c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
    result = qt.mcsolve(
        driven_hamiltonian(ops, CANDIDATES[sector], 0.5), psi0,
        np.linspace(0.0, t_stop, 9), c_ops, ntraj=ntraj, seeds=seed,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    return result.states[-1]


def main():
    out = {}
    for n_max in (2, 3):
        n_cells, sector, profile, disorder = 2, "chiral_topological", "edge_low_loss", 0.0
        draw_seed = 20260814  # sample 0
        state = run_condition(n_cells, n_max, sector, profile, disorder,
                              ntraj=10, seed=draw_seed + 100)
        n_modes = 2 * n_cells
        c = edge_mode_vector(n_cells, CANDIDATES[sector])
        u = complete_mode_rotation(c)
        T = many_body_rotation(u, n_modes, n_max)
        rho_rot = rotate_state(state, T, n_modes, n_max)
        rho_local = state.ptrace(0)
        rho_edge = rho_rot.ptrace(0)
        entry = {"n_max": n_max, "sector": sector, "profile": profile,
                 "disorder": disorder, "draw_seed": draw_seed}
        for gp in (31, 61, 91):
            neg_l, mn_l, mx_l = wigner_stats(rho_local, gp)
            neg_e, mn_e, mx_e = wigner_stats(rho_edge, gp)
            entry[f"grid{gp}"] = {
                "local_negativity": neg_l, "local_wigner_min": mn_l,
                "collective_negativity": neg_e, "collective_wigner_min": mn_e,
            }
        out[f"n_max_{n_max}"] = entry
    payload = {
        "status": "COMPLETED_COLLECTIVE_WIGNER_RESOLUTION_CHECK",
        "note": ("collective edge-mode negativity vs grid resolution and Fock "
                 "truncation; representative topological/edge_low_loss/d=0/sample0"),
        "results": out,
    }
    dest = ROOT / "data" / "collective_wigner_resolution_check.json"
    dest.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"saved {dest}")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
