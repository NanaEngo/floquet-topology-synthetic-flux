#!/usr/bin/env python3
"""Run paired dynamics for the explicit open-chain Floquet mapping hypothesis.

This is intentionally separate from run_matched_dynamics.py, which models the
historical static chain.  It must not overwrite or relabel historical outputs.
The mapped model is a validation pilot, not yet a certified physical realization
of the archived labels.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import qutip as qt

ROOT = Path(__file__).resolve().parents[1]
J_CANDIDATES = {"candidate_W1": (0.2, 1.2), "candidate_W0": (0.1, 0.1)}
MASS = 1.0
DELTA = 0.5
PERIOD = 2.0 * np.pi
GAMMA_EDGE_DEFAULT = 1.0e-4
GAMMA_BULK_DEFAULT = 0.2
ALPHA = 1.0
DISORDER_LEVELS = (0.0, 0.5, 1.0, 1.5)


def build_ops(n_modes: int, n_max: int):
    a = qt.destroy(n_max + 1)
    eye = qt.qeye(n_max + 1)
    return [qt.tensor(*[a if i == site else eye for i in range(n_modes)])
            for site in range(n_modes)]


def single_particle_matrix(n_cells: int, j1: float, j2: float, step: int):
    n = 2 * n_cells
    h = np.zeros((n, n), dtype=complex)
    m_eff = MASS - step * DELTA
    for cell in range(n_cells):
        a, b = 2 * cell, 2 * cell + 1
        h[a, a] = m_eff
        h[b, b] = -m_eff
    for cell in range(n_cells - 1):
        a0, b0 = 2 * cell, 2 * cell + 1
        a1, b1 = 2 * (cell + 1), 2 * (cell + 1) + 1
        for i, j, value in ((a0, a1, j2 / 2.0),
                            (b0, b1, -j2 / 2.0),
                            (a0, b1, -1j * j1 / 2.0),
                            (a1, b0, 1j * j1 / 2.0)):
            h[i, j] = value
            h[j, i] = np.conjugate(value)
    return h


def second_quantize(matrix: np.ndarray, ops):
    h = 0 * ops[0]
    for i in range(matrix.shape[0]):
        if abs(matrix[i, i]) > 0:
            h += float(np.real(matrix[i, i])) * (ops[i].dag() * ops[i])
        for j in range(i + 1, matrix.shape[1]):
            value = matrix[i, j]
            if abs(value) > 0:
                h += value * ops[i].dag() * ops[j]
                h += np.conjugate(value) * ops[j].dag() * ops[i]
    return h


def step_hamiltonian(ops, j1, j2, step):
    h = second_quantize(single_particle_matrix(len(ops) // 2, j1, j2, step), ops)
    h += 0.5 * (ops[0].dag() * ops[0].dag() * ops[0] * ops[0])
    return h


def coeff_step0(t, **kwargs):
    return 1.0 if (float(t) % PERIOD) < PERIOD / 2.0 else 0.0


def coeff_step1(t, **kwargs):
    return 1.0 - coeff_step0(t, **kwargs)


def driven_hamiltonian(ops, j1, j2):
    h0 = step_hamiltonian(ops, j1, j2, 0)
    h1 = step_hamiltonian(ops, j1, j2, 1)
    return qt.QobjEvo([[h0, coeff_step0], [h1, coeff_step1]])


def initial_state(n_modes, n_max):
    dim = n_max + 1
    cat = (qt.coherent(dim, ALPHA) + qt.coherent(dim, -ALPHA)).unit()
    vacuum = qt.basis(dim, 0)
    return qt.tensor(*[cat if i == 0 else vacuum for i in range(n_modes)])


def rates(n_modes, profile, disorder_level, seed, gamma_edge, gamma_bulk):
    base = np.full(n_modes, gamma_bulk, dtype=float)
    if profile == "edge_low_loss":
        base[[0, -1]] = gamma_edge
    elif profile != "uniform_loss":
        raise ValueError(profile)
    rng = np.random.default_rng(seed)
    return np.maximum(base + rng.uniform(-disorder_level, disorder_level) * base, 0.0)


def negativity(state, site, grid_points):
    rho = state.ptrace(site)
    if rho.isket:
        rho = qt.ket2dm(rho)
    grid = np.linspace(-4.0, 4.0, grid_points)
    w = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    return float(np.sum(np.abs(w[w < 0.0])) * dx * dx)


def simulate(h, psi0, c_ops, ops, ntraj, seed, t_stop, n_time, grid_points):
    result = qt.mcsolve(
        h, psi0, np.linspace(0.0, t_stop, n_time), c_ops,
        ntraj=ntraj, seeds=seed,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    state = result.states[-1]
    occ = np.asarray([float(np.real(qt.expect(op.dag() * op, state))) for op in ops])
    return {
        "negativity_edge": negativity(state, 0, grid_points),
        "negativity_center": negativity(state, len(ops) // 2, grid_points),
        "edge_population": float(occ[0] + occ[-1]),
        "bulk_population": float(np.sum(occ[1:-1])),
        "occupations": occ.tolist(),
    }


def run(args):
    n_modes = 2 * args.n_cells
    ops = build_ops(n_modes, args.n_max)
    psi0 = initial_state(n_modes, args.n_max)
    output_rows = []
    start = time.time()
    for profile_index, profile in enumerate(("edge_low_loss", "uniform_loss")):
        for disorder_index, disorder_level in enumerate(DISORDER_LEVELS):
            draw_base = args.seed + profile_index * 1_000_000 + disorder_index * 10_000
            for sample in range(args.n_mc):
                draw_seed = draw_base + sample
                rate_values = rates(n_modes, profile, disorder_level, draw_seed,
                                    args.gamma_edge, args.gamma_bulk)
                c_ops = [np.sqrt(x) * op for x, op in zip(rate_values, ops) if x > 0]
                for sector_index, (sector, (j1, j2)) in enumerate(J_CANDIDATES.items()):
                    row = simulate(
                        driven_hamiltonian(ops, j1, j2), psi0, c_ops, ops,
                        args.ntraj, draw_seed + 100 + sector_index,
                        args.t_stop, args.n_time, args.wigner_grid,
                    )
                    output_rows.append({
                        "n_cells": args.n_cells,
                        "n_modes": n_modes,
                        "n_max": args.n_max,
                        "sector": sector,
                        "provisional_winding": 1 if sector == "candidate_W1" else 0,
                        "J1": j1,
                        "J2": j2,
                        "loss_profile": profile,
                        "disorder_level": disorder_level,
                        "sample": sample,
                        "draw_seed": draw_seed,
                        "trajectory_seed": draw_seed + 100 + sector_index,
                        "ntraj": args.ntraj,
                        "t_stop": args.t_stop,
                        "period": PERIOD,
                        "wigner_grid": args.wigner_grid,
                        "rates": rate_values.tolist(),
                        **row,
                    })
            print(f"completed profile={profile} disorder={disorder_level}")
    output = {
        "status": "COMPLETED_MAPPED_FLOQUET_PILOT",
        "claim_scope": "finite-time dynamics for an explicit open-chain mapping hypothesis",
        "physical_sector_certified": False,
        "mapping_document": "data/model_mapping_hypothesis.md",
        "not_supported": ["topological protection theorem", "certified archived-label mapping", "experimental validation", "metrological gain"],
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "design": {
            "n_cells": args.n_cells, "n_modes": n_modes, "n_max": args.n_max,
            "n_mc": args.n_mc, "ntraj": args.ntraj, "wigner_grid": args.wigner_grid,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": args.t_stop, "n_time_points": args.n_time,
            "period": PERIOD, "paired_disorder": True,
            "independent_unit": "disorder realization; trajectories nested",
            "disorder_levels": DISORDER_LEVELS,
        },
        "records": output_rows,
        "elapsed_seconds": time.time() - start,
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(output, indent=2) + "\n")
    print(f"saved {out} with {len(output_rows)} records")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-cells", type=int, default=2)
    p.add_argument("--n-max", type=int, default=2)
    p.add_argument("--n-mc", type=int, default=2)
    p.add_argument("--ntraj", type=int, default=5)
    p.add_argument("--wigner-grid", type=int, default=31)
    p.add_argument("--gamma-edge", type=float, default=GAMMA_EDGE_DEFAULT)
    p.add_argument("--gamma-bulk", type=float, default=GAMMA_BULK_DEFAULT)
    p.add_argument("--t-stop", type=float, default=PERIOD)
    p.add_argument("--n-time", type=int, default=9)
    p.add_argument("--seed", type=int, default=20260814)
    p.add_argument("--output", type=Path, default=ROOT / "data" / "mapped_floquet_pilot.json")
    args = p.parse_args()
    if min(args.n_cells, args.n_max, args.n_mc, args.ntraj) < 1:
        p.error("counts must be positive")
    if args.t_stop <= 0 or args.n_time < 2 or args.wigner_grid < 11:
        p.error("invalid time or grid")
    run(args)


if __name__ == "__main__":
    main()
