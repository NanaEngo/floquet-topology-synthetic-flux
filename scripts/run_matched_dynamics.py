#!/usr/bin/env python3
"""Run paired finite-time dynamics for provisional source-model sector candidates.

The two candidates are selected only after the FL_QOM resolution/gap audit.  They
remain *source-model candidates*, not certified physical topological sectors,
until the mapping between the FL_QOM invariant model and this open N-site model
is resolved.

The independent unit is one disorder realization.  Quantum trajectories are
nested within that realization and are not treated as independent observations.
The same disorder draw is used for both candidates within a profile.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import qutip as qt

ROOT = Path(__file__).resolve().parents[1]
CANDIDATE_FILE = ROOT / "data" / "candidate_sector_validation.json"
N_MAX_DEFAULT = 2
J_CANDIDATES = {"candidate_W1": (0.2, 1.2), "candidate_W0": (0.1, 0.1)}
GAMMA_EDGE = 1.0e-4
GAMMA_BULK = 1.0
ALPHA = 1.0
DISORDER_LEVELS = (0.0, 0.5, 1.0, 1.5)
T_STOP = 1.0
N_TIME = 6


def build_ops(n_sites: int, n_max: int):
    a = qt.destroy(n_max + 1)
    eye = qt.qeye(n_max + 1)
    ops = []
    for site in range(n_sites):
        ops.append(qt.tensor(*[a if i == site else eye for i in range(n_sites)]))
    return ops


def build_hamiltonian(ops, j1: float, j2: float) -> qt.Qobj:
    h = 0 * ops[0]
    for i in range(len(ops) - 1):
        hopping = j1 if i % 2 == 0 else j2
        h += hopping * (ops[i].dag() * ops[i + 1] + ops[i] * ops[i + 1].dag())
    h += 0.5 * (ops[0].dag() * ops[0].dag() * ops[0] * ops[0])
    return h


def initial_state(n_sites: int, n_max: int) -> qt.Qobj:
    dim = n_max + 1
    cat = (qt.coherent(dim, ALPHA) + qt.coherent(dim, -ALPHA)).unit()
    vacuum = qt.basis(dim, 0)
    return qt.tensor(*[cat if i == 0 else vacuum for i in range(n_sites)])


def rates_for_profile(n_sites: int, profile: str, disorder_level: float, seed: int,
                       gamma_edge: float, gamma_bulk: float):
    if profile == "edge_low_loss":
        base = np.full(n_sites, gamma_bulk, dtype=float)
        base[[0, -1]] = gamma_edge
    elif profile == "uniform_loss":
        base = np.full(n_sites, gamma_bulk, dtype=float)
    else:
        raise ValueError(f"unknown loss profile: {profile}")
    rng = np.random.default_rng(seed)
    perturbation = rng.uniform(-disorder_level, disorder_level) * base
    return np.maximum(base + perturbation, 0.0)


def collapse_ops(ops, rates):
    return [np.sqrt(rate) * op for op, rate in zip(ops, rates) if rate > 0.0]


def negativity(state: qt.Qobj, site: int, grid_points: int) -> float:
    rho = state.ptrace(site)
    if rho.isket:
        rho = qt.ket2dm(rho)
    grid = np.linspace(-4.0, 4.0, grid_points)
    w = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    return float(np.sum(np.abs(w[w < 0.0])) * dx * dx)


def site_occupations(state: qt.Qobj, ops):
    return np.asarray([float(np.real(qt.expect(op.dag() * op, state))) for op in ops])


def simulate(hamiltonian, psi0, c_ops, ops, ntraj, seed, tlist, grid_points):
    result = qt.mcsolve(
        hamiltonian,
        psi0,
        tlist,
        c_ops,
        ntraj=ntraj,
        seeds=seed,
        options={"nsteps": 5000, "store_states": True, "progress_bar": False},
    )
    state = result.states[-1]
    occ = site_occupations(state, ops)
    return {
        "negativity_edge": negativity(state, 0, grid_points),
        "negativity_center": negativity(state, len(ops) // 2, grid_points),
        "edge_population": float(occ[0] + occ[-1]),
        "bulk_population": float(np.sum(occ[1:-1])),
        "occupations": occ.tolist(),
    }


def verify_candidates():
    data = json.loads(CANDIDATE_FILE.read_text())
    for name, expected in (("candidate_W1", 1), ("candidate_W0", 0)):
        rec = data["candidates"][name]
        if not rec["source_model_candidate_pass"]:
            raise RuntimeError(f"{name} failed source-model candidate audit")
        if rec["provisional_label"] != expected:
            raise RuntimeError(f"unexpected provisional label for {name}")


def run(args):
    verify_candidates()
    ops = build_ops(args.n_sites, args.n_max)
    psi0 = initial_state(args.n_sites, args.n_max)
    tlist = np.linspace(0.0, args.t_stop, args.n_time)
    outputs = []
    start = time.time()
    for profile_index, profile in enumerate(("edge_low_loss", "uniform_loss")):
        for disorder_index, disorder_level in enumerate(DISORDER_LEVELS):
            for sample in range(args.n_mc):
                draw_seed = args.seed + profile_index * 1_000_000 + disorder_index * 10_000 + sample
                rates = rates_for_profile(args.n_sites, profile, disorder_level, draw_seed,
                                          args.gamma_edge, args.gamma_bulk)
                c_ops = collapse_ops(ops, rates)
                for sector_index, (sector, (j1, j2)) in enumerate(J_CANDIDATES.items()):
                    run_seed = draw_seed + 100 + sector_index
                    values = simulate(
                        build_hamiltonian(ops, j1, j2), psi0, c_ops, ops,
                        args.ntraj, run_seed, tlist, args.wigner_grid
                    )
                    outputs.append({
                        "n_sites": args.n_sites,
                        "n_max": args.n_max,
                        "sector": sector,
                        "provisional_winding": 1 if sector == "candidate_W1" else 0,
                        "J1": j1,
                        "J2": j2,
                        "loss_profile": profile,
                        "disorder_level": disorder_level,
                        "sample": sample,
                        "draw_seed": draw_seed,
                        "trajectory_seed": run_seed,
                        "ntraj": args.ntraj,
                        "t_stop": args.t_stop,
                        "wigner_grid": args.wigner_grid,
                        "rates": rates.tolist(),
                        **values,
                    })
            print(f"completed profile={profile} disorder={disorder_level}")
    output = {
        "status": "COMPLETED_SOURCE_MODEL_PAIRED_DYNAMICS",
        "claim_scope": "finite-time paired numerical dynamics conditional on audited source-model candidates",
        "physical_sector_certified": False,
        "not_supported": ["topological protection theorem", "certified open-system invariant", "experimental validation", "metrological gain"],
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "design": {
            "n_sites": args.n_sites, "n_max": args.n_max, "n_mc": args.n_mc,
            "ntraj": args.ntraj, "wigner_grid": args.wigner_grid,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": args.t_stop, "n_time_points": args.n_time,
            "disorder_levels": DISORDER_LEVELS,
            "loss_profiles": ["edge_low_loss", "uniform_loss"],
            "independent_unit": "disorder realization; trajectories nested",
            "paired_disorder": True,
            "seed": args.seed,
        },
        "candidate_source": "data/candidate_sector_validation.json",
        "records": outputs,
        "elapsed_seconds": time.time() - start,
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(output, indent=2) + "\n")
    print(f"saved {out} with {len(outputs)} records")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-sites", type=int, default=6)
    parser.add_argument("--n-max", type=int, default=N_MAX_DEFAULT)
    parser.add_argument("--n-mc", type=int, default=10)
    parser.add_argument("--ntraj", type=int, default=10)
    parser.add_argument("--wigner-grid", type=int, default=31)
    parser.add_argument("--gamma-edge", type=float, default=GAMMA_EDGE)
    parser.add_argument("--gamma-bulk", type=float, default=GAMMA_BULK)
    parser.add_argument("--t-stop", type=float, default=T_STOP)
    parser.add_argument("--n-time", type=int, default=N_TIME)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "matched_dynamics_N6_pilot.json")
    args = parser.parse_args()
    if (args.n_sites < 4 or args.n_max < 1 or args.n_mc < 1 or args.ntraj < 1
            or args.wigner_grid < 11 or args.n_time < 2
            or args.gamma_edge < 0 or args.gamma_bulk < 0 or args.t_stop <= 0):
        parser.error("invalid size, sampling, grid, rate or time arguments")
    run(args)


if __name__ == "__main__":
    main()
