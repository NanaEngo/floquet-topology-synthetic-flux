#!/usr/bin/env python3
"""Reviewer-requested convergence checks for the chiral Floquet production.

1. Local-negativity contrast vs Fock cutoff N_max in {2,3,4} (2 cells, disorder 0).
2. Wigner-grid convergence of the single-site negativity integral on a fixed state.
3. Trajectory-count sensitivity of the local negativity estimate (nested within one
   disorder realization; trajectories are not independent units).

All parameters match the production protocol (run_chiral_floquet_dynamics.py):
period 2*pi, kerr 0.5, gamma_edge 1e-4, gamma_bulk 0.2, ntraj 10 (production),
even cat preparation on the boundary site, edge_site initial-state mode.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import qutip as qt
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
PERIOD = 2.0 * np.pi
HALF_PERIOD = np.pi
CANDIDATES = {
    "chiral_topological": {"intra_A": 0.1, "inter_B": 1.2},
    "chiral_trivial": {"intra_A": 1.2, "inter_B": 0.1},
}
ALPHA = 1.0
GAMMA_EDGE = 1e-4
GAMMA_BULK = 0.2
KERR = 0.5
N_CELLS = 2
N_MC = 20
NTRAJ = 10
SEED_BASE = 20260815


def build_ops(n_modes, n_max):
    a = qt.destroy(n_max + 1)
    eye = qt.qeye(n_max + 1)
    return [qt.tensor(*[a if i == j else eye for i in range(n_modes)])
            for j in range(n_modes)]


def finite_chain_matrix(n_cells, intra, inter):
    n = 2 * n_cells
    h = np.zeros((n, n), dtype=complex)
    for cell in range(n_cells):
        a, b = 2 * cell, 2 * cell + 1
        h[a, b] = intra
        h[b, a] = intra
        if cell < n_cells - 1:
            an, bn = 2 * (cell + 1), 2 * (cell + 1) + 1
            h[a, bn] = inter
            h[bn, a] = inter
    return h


def second_quantize(matrix, ops):
    h = 0 * ops[0]
    for i in range(matrix.shape[0]):
        for j in range(i + 1, matrix.shape[1]):
            if abs(matrix[i, j]) > 0:
                h += matrix[i, j] * ops[i].dag() * ops[j]
                h += np.conjugate(matrix[i, j]) * ops[j].dag() * ops[i]
    return h


def step_hamiltonian(ops, intra, inter, kerr):
    h = second_quantize(finite_chain_matrix(len(ops) // 2, intra, inter), ops)
    h += kerr * (ops[0].dag() * ops[0].dag() * ops[0] * ops[0])
    return h


def coeff_A(t, **kwargs):
    return 1.0 if (float(t) % PERIOD) < HALF_PERIOD else 0.0


def coeff_B(t, **kwargs):
    return 1.0 - coeff_A(t, **kwargs)


def driven_hamiltonian(ops, candidate, kerr):
    ha = step_hamiltonian(ops, candidate["intra_A"], 0.0, kerr)
    hb = step_hamiltonian(ops, 0.0, candidate["inter_B"], kerr)
    return qt.QobjEvo([[ha, coeff_A], [hb, coeff_B]])


def initial_edge_site(n_cells, n_max):
    dim = n_max + 1
    n_modes = 2 * n_cells
    amplitudes = np.zeros(n_modes, dtype=complex)
    amplitudes[0] = ALPHA
    plus = qt.tensor(*[qt.coherent(dim, value) for value in amplitudes])
    minus = qt.tensor(*[qt.coherent(dim, -value) for value in amplitudes])
    return (plus + minus).unit()


def negativity(state, grid_points):
    rho = state.ptrace(0)
    if rho.isket:
        rho = qt.ket2dm(rho)
    grid = np.linspace(-4.0, 4.0, grid_points)
    w = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    return float(np.sum(np.abs(w[w < 0])) * dx * dx)


def simulate(h, psi0, c_ops, ntraj, seed, grid_points=31):
    result = qt.mcsolve(
        h, psi0, np.linspace(0.0, PERIOD, 9), c_ops,
        ntraj=ntraj, seeds=seed,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    return negativity(result.states[-1], grid_points)


def run(args, out_path):
    start = time.time()
    rows = []
    for n_max in (2, 3, 4):
        n_modes = 2 * N_CELLS
        ops = build_ops(n_modes, n_max)
        psi0 = initial_edge_site(N_CELLS, n_max)
        for profile in ("edge_low_loss", "uniform_loss"):
            base = np.full(n_modes, GAMMA_BULK, dtype=float)
            if profile == "edge_low_loss":
                base[[0, -1]] = GAMMA_EDGE
            for sample in range(N_MC):
                draw_seed = args.seed + n_max * 1_000_000 + (0 if profile == "edge_low_loss" else 500_000) + sample
                rng = np.random.default_rng(draw_seed)
                rates = np.maximum(base + rng.uniform(-0.0, 0.0) * base, 0.0)  # disorder 0.0
                c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
                for sector, candidate in CANDIDATES.items():
                    seed = draw_seed + 100 + (0 if sector == "chiral_topological" else 1)
                    neg = simulate(driven_hamiltonian(ops, candidate, KERR), psi0, c_ops, NTRAJ, seed)
                    rows.append({
                        "n_max": n_max, "n_modes": n_modes, "loss_profile": profile,
                        "disorder_level": 0.0, "sample": sample, "sector": sector,
                        "negativity_edge": neg, "ntraj": NTRAJ,
                    })
        print(f"completed n_max={n_max}", flush=True)

    # Wigner-grid convergence on a fixed final state (first topological row of n_max=2)
    ref = next(r for r in rows if r["n_max"] == 2 and r["loss_profile"] == "edge_low_loss"
               and r["sector"] == "chiral_topological" and r["sample"] == 0)
    n_modes = 4
    ops2 = build_ops(n_modes, 2)
    psi0 = initial_edge_site(2, 2)
    base = np.full(n_modes, GAMMA_BULK, dtype=float)
    base[[0, -1]] = GAMMA_EDGE
    rng = np.random.default_rng(SEED_BASE + 7)
    rates = base
    c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops2) if x > 0]
    h = driven_hamiltonian(ops2, CANDIDATES["chiral_topological"], KERR)
    result = qt.mcsolve(h, psi0, np.linspace(0.0, PERIOD, 9), c_ops, ntraj=200, seeds=SEED_BASE + 8,
                        options={"nsteps": 10000, "store_states": True, "progress_bar": False})
    state = result.states[-1]
    rho = state.ptrace(0)
    if rho.isket:
        rho = qt.ket2dm(rho)
    grid_rows = []
    for gp in (21, 31, 41, 51, 61):
        grid = np.linspace(-4.0, 4.0, gp)
        w = qt.wigner(rho, grid, grid)
        dx = grid[1] - grid[0]
        grid_rows.append({"grid_points": gp, "negativity_edge": float(np.sum(np.abs(w[w < 0])) * dx * dx)})
    grid_ci = grid_rows[-1]["negativity_edge"]  # finest reference

    # Trajectory-count sensitivity: same realization, ntraj in {10, 50, 200}
    traj_rows = []
    for ntraj in (10, 50, 200):
        val = simulate(h, psi0, c_ops, ntraj, SEED_BASE + 8)
        traj_rows.append({"ntraj": ntraj, "negativity_edge": val})

    output = {
        "status": "CONVERGENCE_CHECKS_COMPLETED",
        "claim_scope": "reviewer-requested descriptive convergence diagnostics only",
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "design": {
            "n_cells": N_CELLS, "kerr": KERR, "gamma_edge": GAMMA_EDGE,
            "gamma_bulk": GAMMA_BULK, "period": PERIOD, "n_mc": N_MC,
            "ntraj_production": NTRAJ, "wigner_grid_production": 31,
            "disorder_levels_tested": [0.0],
            "independent_unit": "disorder realization; trajectories nested",
        },
        "records": rows,
        "grid_convergence": grid_rows,
        "grid_reference_negativity": grid_ci,
        "trajectory_count_sensitivity": traj_rows,
        "elapsed_seconds": time.time() - start,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(output, indent=2) + "\n")
    print(f"saved {out_path} with {len(rows)} records in {output['elapsed_seconds']:.0f} s", flush=True)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=SEED_BASE)
    p.add_argument("--output", type=Path, default=ROOT / "data" / "convergence_checks.json")
    a = p.parse_args()
    run(a, a.output)
