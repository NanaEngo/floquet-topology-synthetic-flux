#!/usr/bin/env python3
"""Matched-design parallel recorder for audit runs R3 (size 5) and R8 (witness).

Extends the production protocol (``run_chiral_floquet_dynamics.py``) without
modifying it.  Every record reproduces the production design bit-for-bit for
the same (profile, disorder, sample, sector): production draw/trajectory seed
scheme, ``n_time=9`` Wigner times, 31x31 grid, ``ntraj=10``, edge_site even-cat
preparation, Kerr chi=0.5.

Run R3 (AUDIT_RESEARCH_REFINEMENTS_20261001.md line 121): fifth size,
5 cells / 10 modes, N_MC=300, disorders {0, 0.5}, both loss profiles, both
sectors.  The size-5 mean is scored against the prediction frozen in
``data/r3_size5_prediction.json`` (L3 out-of-sample test).

Run R8 (line 126): bit-exact rerun of the production design at 2/3/4 cells
with ``--dump-blocks``, storing the terminal edge and centre reduced density
matrices (3x3) for post-hoc certification (Chabaud Fock witnesses and
relative-entropy non-Gaussianity, computed by ``analyze_r8_witness.py``).
The runner also recomputes the production 31x31 negativities so the analyzer
can fail closed if the rerun diverges from the archived production records.

Args (CLI):
    --n-cells           chain length in unit cells (dimensionless)
    --n-mc              number of disorder draws (dimensionless)
    --disorder-levels   multiplicative disorder amplitudes (dimensionless)
    --dump-blocks       also store terminal 3x3 reduced blocks (R8 flag)
    --n-workers         Pool size (dimensionless); pin OMP/BLAS threads to 1
    --output            destination JSON path

Output:
    JSON with production-identical record fields plus optional
    ``rho_edge`` / ``rho_center`` blocks (lists of [Re, Im]).
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
    negativity,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SEED = 20260814  # production seed base


def reduced_block(state: qt.Qobj, site: int) -> list:
    """Extract the terminal single-site reduced density matrix as JSON pairs.

    Args:
        state (qutip.Qobj): Full $n$-mode terminal state (ket or dm).
        site (int): Site index (dimensionless).

    Returns:
        list: $3\times 3$ matrix as rows of ``[Re, Im]`` pairs (dimensionless).
    """
    rho = state.ptrace(site)
    if rho.isket:
        rho = qt.ket2dm(rho)
    return [[[float(x.real), float(x.imag)] for x in row] for row in rho.full()]


def simulate_single(task: tuple) -> dict:
    """Run one (profile, disorder, sample, sector) record. Pool worker.

    Args:
        task (tuple): ``(profile, disorder, sample, sector_key, n_cells, n_max,
            ntraj, wigner_grid, gamma_edge, gamma_bulk, t_stop, n_time, kerr,
            base_seed, dump_blocks)``.

    Returns:
        dict: Production-identical record row, plus ``rho_edge`` and
            ``rho_center`` terminal blocks when ``dump_blocks`` is set.
    """
    (profile, disorder, sample, sector_key, n_cells, n_max, ntraj,
     wigner_grid, gamma_edge, gamma_bulk, t_stop, n_time, kerr,
     base_seed, dump_blocks) = task

    n_modes = 2 * n_cells
    ops = build_ops(n_modes, n_max)  # rebuilt in worker: Qobj not picklable
    candidate = CANDIDATES[sector_key]
    psi0 = initial_state(n_cells, n_max, "edge_site")

    # production index basis keeps draw seeds identical to the archived runs
    disorder_index = DISORDER_LEVELS.index(disorder)
    profile_index = ("edge_low_loss", "uniform_loss").index(profile)
    draw_seed = base_seed + profile_index * 1_000_000 + disorder_index * 10_000 + sample
    rates = rate_vector(n_modes, profile, disorder, draw_seed, gamma_edge, gamma_bulk)
    c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]

    sector_index = list(CANDIDATES.keys()).index(sector_key)
    seed = draw_seed + 100 + sector_index

    h = driven_hamiltonian(ops, candidate, kerr)
    tlist = np.linspace(0.0, t_stop, n_time)
    result = qt.mcsolve(
        h, psi0, tlist, c_ops,
        ntraj=ntraj, seeds=seed,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    state = result.states[-1]
    occ = np.asarray([float(np.real(qt.expect(op.dag() * op, state))) for op in ops])

    row = {
        "n_cells": n_cells, "n_modes": n_modes, "n_max": n_max,
        "sector": sector_key, **candidate, "kerr": kerr,
        "loss_profile": profile, "disorder_level": disorder,
        "sample": sample, "draw_seed": draw_seed,
        "trajectory_seed": seed, "ntraj": ntraj,
        "t_stop": t_stop, "period": PERIOD,
        "wigner_grid": wigner_grid, "rates": rates.tolist(),
        "negativity_edge": negativity(state, 0, wigner_grid),
        "negativity_center": negativity(state, n_modes // 2, wigner_grid),
        "edge_population": float(occ[0] + occ[-1]),
        "bulk_population": float(np.sum(occ[1:-1])),
        "occupations": occ.tolist(),
    }
    if dump_blocks:
        row["rho_edge"] = reduced_block(state, 0)
        row["rho_center"] = reduced_block(state, n_modes // 2)
    return row


def main() -> None:
    """Parse CLI, dispatch tasks over a process pool, write the JSON."""
    p = argparse.ArgumentParser(description="R3/R8 matched-design parallel recorder")
    p.add_argument("--n-cells", type=int, required=True)
    p.add_argument("--n-mc", type=int, default=300)
    p.add_argument("--n-max", type=int, default=2)
    p.add_argument("--disorder-levels", type=float, nargs="+",
                   default=list(DISORDER_LEVELS))
    p.add_argument("--ntraj", type=int, default=10)
    p.add_argument("--wigner-grid", type=int, default=31)
    p.add_argument("--kerr", type=float, default=0.5)
    p.add_argument("--gamma-edge", type=float, default=1e-4)
    p.add_argument("--gamma-bulk", type=float, default=0.2)
    p.add_argument("--t-stop", type=float, default=PERIOD)
    p.add_argument("--n-time", type=int, default=9)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--n-workers", type=int, default=16)
    p.add_argument("--dump-blocks", action="store_true")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    levels = tuple(float(d) for d in args.disorder_levels)
    for d in levels:
        if d not in DISORDER_LEVELS:
            p.error(f"disorder level {d} not in production basis {DISORDER_LEVELS}")

    n_tasks = 2 * len(levels) * args.n_mc * len(CANDIDATES)
    print(f"tasks={n_tasks} workers={args.n_workers} cells={args.n_cells} "
          f"n_mc={args.n_mc} disorders={levels} dump_blocks={args.dump_blocks}",
          flush=True)

    tasks = []
    for profile in ("edge_low_loss", "uniform_loss"):
        for disorder in levels:
            for sample in range(args.n_mc):
                for sector_key in CANDIDATES:
                    tasks.append((
                        profile, disorder, sample, sector_key,
                        args.n_cells, args.n_max, args.ntraj, args.wigner_grid,
                        args.gamma_edge, args.gamma_bulk, args.t_stop,
                        args.n_time, args.kerr, args.seed, args.dump_blocks,
                    ))

    start = time.time()
    rows = []
    with Pool(processes=args.n_workers) as pool:
        for i, row in enumerate(pool.imap_unordered(simulate_single, tasks)):
            rows.append(row)
            if (i + 1) % 200 == 0:
                elapsed = time.time() - start
                rate = (i + 1) / elapsed
                eta = (n_tasks - i - 1) / rate / 60
                print(f"  [{i+1}/{n_tasks}] {rate:.2f} tasks/s  ETA {eta:.0f} min",
                      flush=True)

    elapsed = time.time() - start
    print(f"Completed {len(rows)} tasks in {elapsed:.0f}s "
          f"({len(rows) / elapsed:.2f} tasks/s)", flush=True)

    output = {
        "status": "COMPLETED_R3_R8_RECORDS",
        "audit_source": "AUDIT_RESEARCH_REFINEMENTS_20261001.md",
        "claim_scope": "finite-time descriptive matched dynamics only; no protection claim",
        "not_supported": ["topological protection theorem", "experimental validation",
                          "metrological gain"],
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "design": {
            "n_cells": args.n_cells, "n_modes": 2 * args.n_cells,
            "n_max": args.n_max, "initial_state_mode": "edge_site",
            "kerr": args.kerr, "n_mc": args.n_mc, "ntraj": args.ntraj,
            "wigner_grid": args.wigner_grid,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": args.t_stop, "n_time_points": args.n_time,
            "period": PERIOD, "disorder_levels": list(levels),
            "paired_disorder": True,
            "independent_unit": "disorder realization; trajectories nested",
            "seed": args.seed, "dump_blocks": args.dump_blocks,
            "production_seed_parity": args.seed == DEFAULT_SEED,
        },
        "records": rows,
        "elapsed_seconds": elapsed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(f"saved {args.output} with {len(rows)} records", flush=True)


if __name__ == "__main__":
    main()
