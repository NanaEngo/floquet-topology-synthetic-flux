#!/usr/bin/env python3
"""Long-time stroboscopic contrast recorder for audit run R9.

Answers the referee question behind R9 (AUDIT_RESEARCH_REFINEMENTS_20261001.md
line 127): is the production sign change of the edge negativity contrast at
$t=T$ an artifact of the stroboscopic snapshot, or does it persist on the way
to the stationary state?

Records are produced by the canonical production propagator
``simulate(series=True)`` of ``run_chiral_floquet_dynamics.py`` (Monte Carlo
wave function, ``ntraj=10``, trajectories nested), so the solver, seeds and
record fields are identical to production; only the horizon is extended,
$t_{\rm stop} = H T$ with ``n_time = H + 1``, which makes the stored grid
strictly stroboscopic ($t_k = k T$, spacing exactly one period) and brackets
the production snapshot $t=T$.  The exact-Lindblad alternative was profiled
and rejected: for 3 cells ``mesolve`` propagates a $729^2$ density matrix and
costs 156 s per 8 periods versus 0.75 s for 2 cells; the Monte Carlo ket
propagation is linear in the horizon and stays within the R9 budget.  Solver
correctness is cross-validated in ``data/cross_solver_benchmark_3cells.json``
(fidelity 0.9996).

Draw and trajectory seeds follow the production scheme (base seed 20260814),
so rates and trajectory ensembles of record ``k`` equal the corresponding
production draw.

Args (CLI):
    --n-cells           chain length in unit cells (dimensionless)
    --n-mc              number of disorder draws (dimensionless)
    --horizon-periods   stroboscopic horizon $H$ (dimensionless)
    --ntraj             Monte Carlo trajectories per record (dimensionless)
    --n-workers         Pool size (dimensionless); pin OMP/BLAS threads to 1
    --output            destination JSON path

Output:
    JSON with per-record ``time_series`` of edge/centre Wigner negativity and
    edge/bulk populations on the stroboscopic grid (dimensionless).
"""
from __future__ import annotations

import argparse
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import qutip as qt

from run_chiral_floquet_dynamics import (
    CANDIDATES,
    DISORDER_LEVELS,
    PERIOD,
    build_ops,
    driven_hamiltonian,
    initial_state,
    rate_vector,
    simulate,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SEED = 20260814  # production seed base


def simulate_single(task: tuple) -> dict:
    """Propagate one (profile, disorder, sample, sector) record to the horizon.

    Args:
        task (tuple): ``(profile, disorder, sample, sector_key, n_cells, n_max,
            ntraj, wigner_grid, gamma_edge, gamma_bulk, horizon, kerr, base_seed)``.

    Returns:
        dict: Production-identical record plus the ``time_series`` block
            sampled on the stroboscopic grid $t_k = k T$.
    """
    (profile, disorder, sample, sector_key, n_cells, n_max, ntraj,
     wigner_grid, gamma_edge, gamma_bulk, horizon, kerr, base_seed) = task

    n_modes = 2 * n_cells
    ops = build_ops(n_modes, n_max)
    candidate = CANDIDATES[sector_key]
    psi0 = initial_state(n_cells, n_max, "edge_site")

    disorder_index = DISORDER_LEVELS.index(disorder)
    profile_index = ("edge_low_loss", "uniform_loss").index(profile)
    draw_seed = base_seed + profile_index * 1_000_000 + disorder_index * 10_000 + sample
    rates = rate_vector(n_modes, profile, disorder, draw_seed, gamma_edge, gamma_bulk)
    c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]

    sector_index = list(CANDIDATES.keys()).index(sector_key)
    seed = draw_seed + 100 + sector_index

    t_stop = horizon * PERIOD
    n_time = horizon + 1  # linspace(0, H*T, H+1) has spacing exactly T
    row = simulate(
        driven_hamiltonian(ops, candidate, kerr), psi0, c_ops, ops,
        ntraj, seed, t_stop, n_time, wigner_grid, series=True,
    )
    return {
        "n_cells": n_cells, "n_modes": n_modes, "n_max": n_max,
        "sector": sector_key, **candidate, "kerr": kerr,
        "loss_profile": profile, "disorder_level": disorder,
        "sample": sample, "draw_seed": draw_seed,
        "trajectory_seed": seed, "ntraj": ntraj,
        "t_stop": t_stop, "period": PERIOD, "horizon_periods": horizon,
        "wigner_grid": wigner_grid, "rates": rates.tolist(),
        **row,
    }


def main() -> None:
    """Parse CLI, dispatch tasks over a process pool, write the JSON."""
    p = argparse.ArgumentParser(description="R9 long-time stroboscopic recorder")
    p.add_argument("--n-cells", type=int, required=True)
    p.add_argument("--n-mc", type=int, default=100)
    p.add_argument("--n-max", type=int, default=2)
    p.add_argument("--disorder-levels", type=float, nargs="+",
                   default=list(DISORDER_LEVELS))
    p.add_argument("--horizon-periods", type=int, default=64)
    p.add_argument("--ntraj", type=int, default=10)
    p.add_argument("--wigner-grid", type=int, default=31)
    p.add_argument("--kerr", type=float, default=0.5)
    p.add_argument("--gamma-edge", type=float, default=1e-4)
    p.add_argument("--gamma-bulk", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--n-workers", type=int, default=16)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    levels = tuple(float(d) for d in args.disorder_levels)
    for d in levels:
        if d not in DISORDER_LEVELS:
            p.error(f"disorder level {d} not in production basis {DISORDER_LEVELS}")

    n_tasks = 2 * len(levels) * args.n_mc * len(CANDIDATES)
    print(f"tasks={n_tasks} workers={args.n_workers} cells={args.n_cells} "
          f"n_mc={args.n_mc} horizon={args.horizon_periods}T disorders={levels}",
          flush=True)

    tasks = []
    for profile in ("edge_low_loss", "uniform_loss"):
        for disorder in levels:
            for sample in range(args.n_mc):
                for sector_key in CANDIDATES:
                    tasks.append((
                        profile, disorder, sample, sector_key,
                        args.n_cells, args.n_max, args.ntraj, args.wigner_grid,
                        args.gamma_edge, args.gamma_bulk,
                        args.horizon_periods, args.kerr, args.seed,
                    ))

    start = time.time()
    rows = []
    with Pool(processes=args.n_workers) as pool:
        for i, row in enumerate(pool.imap_unordered(simulate_single, tasks)):
            rows.append(row)
            if (i + 1) % 100 == 0:
                elapsed = time.time() - start
                rate = (i + 1) / elapsed
                eta = (n_tasks - i - 1) / rate / 60
                print(f"  [{i+1}/{n_tasks}] {rate:.2f} tasks/s  ETA {eta:.0f} min",
                      flush=True)

    elapsed = time.time() - start
    print(f"Completed {len(rows)} tasks in {elapsed:.0f}s "
          f"({len(rows) / elapsed:.2f} tasks/s)", flush=True)

    output = {
        "status": "COMPLETED_R9_LONGTIME",
        "audit_source": "AUDIT_RESEARCH_REFINEMENTS_20261001.md",
        "claim_scope": "finite-time descriptive matched dynamics only; no protection claim",
        "not_supported": ["topological protection theorem", "experimental validation",
                          "metrological gain"],
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "design": {
            "n_cells": args.n_cells, "n_modes": 2 * args.n_cells,
            "n_max": args.n_max, "initial_state_mode": "edge_site",
            "kerr": args.kerr, "n_mc": args.n_mc, "ntraj": args.ntraj,
            "solver": "mcsolve (production parity; cross-validated vs mesolve)",
            "wigner_grid": args.wigner_grid,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "horizon_periods": args.horizon_periods, "period": PERIOD,
            "time_grid": "stroboscopic t_k = k*PERIOD (n_time = H+1)",
            "disorder_levels": list(levels), "paired_disorder": True,
            "independent_unit": "disorder realization; trajectories nested",
            "seed": args.seed, "production_seed_parity": args.seed == DEFAULT_SEED,
        },
        "records": rows,
        "elapsed_seconds": elapsed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(f"saved {args.output} with {len(rows)} records", flush=True)


if __name__ == "__main__":
    main()
