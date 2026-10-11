#!/usr/bin/env python3
"""Cross-solver benchmark: QuTiP ``mcsolve`` (Monte Carlo) vs ``mesolve`` (master equation).

Roadmap gate G5, item 6: an independent-solver cross-check on a bounded subset
of the chiral Floquet SSH model. Uses the production parameters at 2 cells /
4 modes with ``n_max=2``, one disorder draw per condition, both certified
sectors and both loss profiles. ``mesolve`` integrates the full Lindblad master
equation (exact in the chosen cutoff); ``mcsolve`` averages ``ntraj`` quantum
trajectories with the same seeds. The final-time averaged states are compared
through Uhlmann fidelity and trace distance, together with the edge population
and the local Wigner negativity used in the main text.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import qutip as qt
from scipy.linalg import expm

from run_chiral_floquet_dynamics import (
    ALPHA,
    CANDIDATES,
    DISORDER_LEVELS,
    HALF_PERIOD,
    PERIOD,
    build_ops,
    driven_hamiltonian,
    finite_chain_matrix,
    initial_state,
    negativity,
    rate_vector,
)

ROOT = Path(__file__).resolve().parents[1]


def trace_distance(rho1, rho2):
    """Trace distance between two density matrices, (1/2)||rho1 - rho2||_1."""
    return float(0.5 * (rho1 - rho2).norm())


def uhlmann_fidelity(rho1, rho2):
    """Uhlmann fidelity F(rho1, rho2), via QuTiP's implementation."""
    if rho1.isket:
        rho1 = qt.ket2dm(rho1)
    if rho2.isket:
        rho2 = qt.ket2dm(rho2)
    return float(qt.fidelity(rho1, rho2))


def run_condition(candidate, n_cells, n_max, kerr, profile, disorder, draw_seed,
                  ntraj, t_stop, n_time, grid_points, gamma_edge, gamma_bulk,
                  trajectory_seed_offset=100):
    n_modes = 2 * n_cells
    ops = build_ops(n_modes, n_max)
    psi0 = initial_state(n_cells, n_max, "edge_site")
    rates = rate_vector(n_modes, profile, disorder, draw_seed, gamma_edge, gamma_bulk)
    c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
    h = driven_hamiltonian(ops, candidate, kerr)
    times = np.linspace(0.0, t_stop, n_time)

    # Exact master equation.
    me = qt.mesolve(
        h, psi0, times, c_ops, e_ops=[],
        options={"nsteps": 10000, "store_states": True},
    )
    rho_me = me.states[-1]

    # Monte Carlo with a fixed trajectory seed.
    mc = qt.mcsolve(
        h, psi0, times, c_ops, e_ops=[],
        ntraj=ntraj, seeds=draw_seed + trajectory_seed_offset,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    rho_mc = mc.states[-1]

    occ_me = np.asarray([float(np.real(qt.expect(op.dag() * op, rho_me))) for op in ops])
    occ_mc = np.asarray([float(np.real(qt.expect(op.dag() * op, rho_mc))) for op in ops])

    return {
        "sector": next(name for name, c in CANDIDATES.items() if c is candidate),
        "n_cells": n_cells, "n_modes": n_modes, "n_max": n_max,
        "kerr": kerr, "loss_profile": profile, "disorder_level": disorder,
        "draw_seed": draw_seed, "ntraj": ntraj,
        "trace_distance": trace_distance(rho_mc, rho_me),
        "uhlmann_fidelity": uhlmann_fidelity(rho_mc, rho_me),
        "edge_population_mc": float(occ_mc[0] + occ_mc[-1]),
        "edge_population_me": float(occ_me[0] + occ_me[-1]),
        "negativity_edge_mc": negativity(rho_mc, 0, grid_points),
        "negativity_edge_me": negativity(rho_me, 0, grid_points),
        "occupation_rmse": float(np.sqrt(np.mean((occ_mc - occ_me) ** 2))),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-cells", type=int, default=2)
    p.add_argument("--n-max", type=int, default=2)
    p.add_argument("--ntraj", type=int, default=200)
    p.add_argument("--wigner-grid", type=int, default=31)
    p.add_argument("--kerr", type=float, default=0.5)
    p.add_argument("--gamma-edge", type=float, default=1e-4)
    p.add_argument("--gamma-bulk", type=float, default=0.2)
    p.add_argument("--t-stop", type=float, default=PERIOD)
    p.add_argument("--n-time", type=int, default=9)
    p.add_argument("--seed", type=int, default=20260814)
    p.add_argument("--disorders", type=float, nargs="+", default=[0.0, 0.5])
    p.add_argument("--output", type=Path, default=ROOT / "data" / "cross_solver_benchmark.json")
    args = p.parse_args()

    rows = []
    start = time.time()
    for profile_index, profile in enumerate(("edge_low_loss", "uniform_loss")):
        for disorder_index, disorder in enumerate(args.disorders):
            draw_seed = args.seed + profile_index * 1_000_000 + disorder_index * 10_000
            for candidate in CANDIDATES.values():
                row = run_condition(
                    candidate, args.n_cells, args.n_max, args.kerr, profile, disorder,
                    draw_seed, args.ntraj, args.t_stop, args.n_time, args.wigner_grid,
                    args.gamma_edge, args.gamma_bulk,
                )
                rows.append(row)
                print(f"done {row['sector']} {profile} d={disorder} "
                      f"F={row['uhlmann_fidelity']:.6f} T={row['trace_distance']:.2e}")

    fidelities = [r["uhlmann_fidelity"] for r in rows]
    traces = [r["trace_distance"] for r in rows]
    rmse = [r["occupation_rmse"] for r in rows]
    output = {
        "status": "COMPLETED_CROSS_SOLVER_BENCHMARK",
        "claim_scope": "independent-solver cross-check on a bounded chiral Floquet SSH subset",
        "design": {
            "n_cells": args.n_cells, "n_modes": 2 * args.n_cells, "n_max": args.n_max,
            "ntraj": args.ntraj, "wigner_grid": args.wigner_grid, "kerr": args.kerr,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": args.t_stop, "n_time_points": args.n_time,
            "period": PERIOD, "disorder_levels": args.disorders,
            "primary_solver": "QuTiP mesolve (exact Lindblad)",
            "cross_solver": "QuTiP mcsolve (Monte Carlo wave function)",
            "comparison_metrics": ["uhlmann_fidelity", "trace_distance", "occupation_rmse"],
            "independent_unit": "disorder draw; trajectories nested",
        },
        "summary": {
            "n_conditions": len(rows),
            "fidelity_min": float(np.min(fidelities)),
            "fidelity_mean": float(np.mean(fidelities)),
            "trace_distance_max": float(np.max(traces)),
            "trace_distance_mean": float(np.mean(traces)),
            "occupation_rmse_max": float(np.max(rmse)),
            "occupation_rmse_mean": float(np.mean(rmse)),
        },
        "records": rows,
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "elapsed_seconds": time.time() - start,
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(output, indent=2) + "\n")
    print(f"saved {out} with {len(rows)} records; "
          f"mean F={output['summary']['fidelity_mean']:.5f} "
          f"max T={output['summary']['trace_distance_max']:.2e}")


if __name__ == "__main__":
    main()
