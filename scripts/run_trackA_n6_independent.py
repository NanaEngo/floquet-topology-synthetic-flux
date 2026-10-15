#!/usr/bin/env python3
"""Independent deterministic N=6 Track-A propagation.

This reproduces the model family used by the archived N=6 bootstrap script while
making the disorder and Monte-Carlo trajectories reproducible.  It writes only to
the V2608 workspace and does not modify the parent repository cache.

The output is a finite-time numerical validation, not evidence for a protection
theorem or an asymptotic bound.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import qutip as qt


N_SITES = 6
N_MAX = 2
DIM = N_MAX + 1
GAMMA_EDGE = 1.0e-4
GAMMA_BULK = 1.0
J_TOPO = (0.1, 1.0)
J_TRIV = (1.0, 1.0)
ALPHA = 1.0
T_STOP = 1.0
N_TIME = 6
DISORDER_LEVELS = (0.1, 0.575, 1.05, 1.525, 2.0)


def operators() -> list[qt.Qobj]:
    a = qt.destroy(DIM)
    eye = qt.qeye(DIM)
    return [qt.tensor(*[a if i == site else eye for i in range(N_SITES)])
            for site in range(N_SITES)]


A_OPS = operators()


def build_hamiltonian(j1: float, j2: float) -> qt.Qobj:
    h = 0 * A_OPS[0]
    for i in range(N_SITES - 1):
        hopping = j1 if i % 2 == 0 else j2
        h += hopping * (A_OPS[i].dag() * A_OPS[i + 1]
                        + A_OPS[i] * A_OPS[i + 1].dag())
    # Match the historical source model's local Kerr term.
    h += 0.5 * (A_OPS[0].dag() * A_OPS[0].dag()
                * A_OPS[0] * A_OPS[0])
    return h


def initial_state() -> qt.Qobj:
    cat = (qt.coherent(DIM, ALPHA) + qt.coherent(DIM, -ALPHA)).unit()
    vacuum = qt.basis(DIM, 0)
    return qt.tensor(*[cat if i == 0 else vacuum for i in range(N_SITES)])


def collapse_operators(disorder: np.ndarray) -> list[qt.Qobj]:
    rates = np.full(N_SITES, GAMMA_BULK, dtype=float)
    rates[[0, -1]] = GAMMA_EDGE
    rates = np.maximum(rates + disorder, 0.0)
    return [np.sqrt(rate) * A_OPS[i] for i, rate in enumerate(rates)]


def wigner_negativity(state: qt.Qobj) -> float:
    rho = state.ptrace(0)
    if rho.isket:
        rho = qt.ket2dm(rho)
    grid = np.linspace(-4.0, 4.0, 31)
    wigner = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    return float(np.sum(np.abs(wigner[wigner < 0.0])) * dx * dx)


def propagate(j1: float, j2: float, collapse: list[qt.Qobj],
              ntraj: int, seed: int, tlist: np.ndarray) -> np.ndarray:
    result = qt.mcsolve(
        build_hamiltonian(j1, j2),
        initial_state(),
        tlist,
        collapse,
        ntraj=ntraj,
        seeds=seed,
        options={"nsteps": 5000, "store_states": True, "progress_bar": False},
    )
    return np.asarray([wigner_negativity(state) for state in result.states], dtype=float)


def run(args: argparse.Namespace) -> dict:
    start = time.time()
    tlist = np.linspace(0.0, T_STOP, N_TIME)
    clean_collapse = collapse_operators(np.zeros(N_SITES))
    clean_topo = propagate(*J_TOPO, clean_collapse, args.clean_ntraj,
                           args.seed, tlist)
    clean_triv = propagate(*J_TRIV, clean_collapse, args.clean_ntraj,
                           args.seed + 1, tlist)

    def decay_rate(values: np.ndarray) -> float:
        mask = values > 1e-15
        if mask.sum() < 2:
            return 0.0
        return float(abs(np.polyfit(tlist[mask], np.log(values[mask]), 1)[0]))

    records: list[dict] = []
    for level_index, disorder_level in enumerate(DISORDER_LEVELS):
        topo_values = []
        triv_values = []
        for sample in range(args.n_mc):
            seed = args.seed + 100_000 * (level_index + 1) + sample
            rng = np.random.default_rng(seed)
            disorder = np.array([
                rng.uniform(-disorder_level, disorder_level)
                * (GAMMA_EDGE if i in (0, N_SITES - 1) else GAMMA_BULK)
                for i in range(N_SITES)
            ])
            collapse = collapse_operators(disorder)
            topo_values.append(propagate(*J_TOPO, collapse, args.ntraj,
                                         seed, tlist)[-1])
            triv_values.append(propagate(*J_TRIV, collapse, args.ntraj,
                                         seed + 1, tlist)[-1])

        topo = np.asarray(topo_values, dtype=float)
        triv = np.asarray(triv_values, dtype=float)
        records.append({
            "disorder_level": float(disorder_level),
            "n_samples": int(args.n_mc),
            "topological_mean": float(topo.mean()),
            "topological_std": float(topo.std(ddof=1)) if len(topo) > 1 else 0.0,
            "trivial_mean": float(triv.mean()),
            "trivial_std": float(triv.std(ddof=1)) if len(triv) > 1 else 0.0,
            "topological_values": topo.tolist(),
            "trivial_values": triv.tolist(),
        })
        print(
            f"D={disorder_level:.3f}: topo={topo.mean():.6f} +/- {topo.std():.6f}; "
            f"trivial={triv.mean():.6f} +/- {triv.std():.6f}"
        )

    output = {
        "status": "COMPLETED",
        "claim_scope": "finite-time deterministic N=6 numerical propagation",
        "not_supported": ["topological protection theorem", "asymptotic bound",
                          "experimental validation", "metrological gain"],
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "model": {
            "n_sites": N_SITES, "n_max": N_MAX, "gamma_edge": GAMMA_EDGE,
            "gamma_bulk": GAMMA_BULK, "J_topological": J_TOPO,
            "J_trivial": J_TRIV, "alpha": ALPHA,
        },
        "sampling": {
            "seed": args.seed, "n_mc": args.n_mc, "ntraj": args.ntraj,
            "clean_ntraj": args.clean_ntraj, "t_stop": T_STOP,
            "n_time_points": N_TIME, "disorder_levels": DISORDER_LEVELS,
        },
        "clean": {
            "topological_trace": clean_topo.tolist(),
            "trivial_trace": clean_triv.tolist(),
            "topological_decay_rate": decay_rate(clean_topo),
            "trivial_decay_rate": decay_rate(clean_triv),
        },
        "disorder_records": records,
        "elapsed_seconds": time.time() - start,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2) + "\n")
    print(f"Saved {output_path}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-mc", type=int, default=500)
    parser.add_argument("--ntraj", type=int, default=10)
    parser.add_argument("--clean-ntraj", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1]
                        / "data" / "trackA_n6_independent.json")
    args = parser.parse_args()
    if min(args.n_mc, args.ntraj, args.clean_ntraj) < 1:
        parser.error("trajectory and sample counts must be positive")
    run(args)


if __name__ == "__main__":
    main()
