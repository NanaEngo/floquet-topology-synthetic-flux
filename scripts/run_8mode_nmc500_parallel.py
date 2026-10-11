#!/usr/bin/env python3
"""Parallelized 8-mode NMC=500 runner for chiral Floquet dynamics.

Splits 8000 independent simulations across CPU cores using multiprocessing.
Each simulation: 1 (profile, disorder, sample, sector) combination.

Physics: the initial state is an even cat ((|+alpha> + |-alpha>)/sqrt(2)),
site occupations then carry only EVEN total quanta; the drive is
number-conserving, the Kerr term is diagonal, and loss only lowers quanta.
The state therefore never leaves the total-quanta Q<=q_max sector, where
q_max is the largest even number <= n_max (Q<=2 in production, Q<=4 for
the n_max=4 cutoff check; see q_max_for). All dynamics are run exactly in
that subspace (dim 45 of (n_max+1)^n_modes, e.g. 45 of
65536 for n_max=3, 8 modes) and observables are embedded back for the Wigner
negativity. Validated 2026-10-02 against the full-space archived rows
(operator-level invariance block diff 0; record-level diffs <= 3e-5, i.e.
integrator branch noise far below the production standard error).

Usage:
    mamba run -n qom python scripts/run_8mode_nmc500_parallel.py \
        --n-cells 4 --n-mc 500 --n-workers 16 \
        --output data/chiral_production/chiral_floquet_4cells_NMC500.json
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import qutip as qt
import scipy.sparse as sp

# Reuse all physics functions from the original script
ROOT = Path(__file__).resolve().parents[1]
PERIOD = 2.0 * np.pi
HALF_PERIOD = np.pi
CANDIDATES = {
    "chiral_topological": {"intra_A": 0.1, "inter_B": 1.2, "nu_0": -1, "nu_pi": 0},
    "chiral_trivial": {"intra_A": 1.2, "inter_B": 0.1, "nu_0": 0, "nu_pi": 0},
}
DISORDER_LEVELS = (0.0, 0.5, 1.0, 1.5)
ALPHA = 1.0


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


def edge_mode_vector(n_cells, candidate):
    from scipy.linalg import expm
    ha = finite_chain_matrix(n_cells, candidate["intra_A"], 0.0)
    hb = finite_chain_matrix(n_cells, 0.0, candidate["inter_B"])
    u = expm(-1j * hb * HALF_PERIOD) @ expm(-1j * ha * HALF_PERIOD)
    values, vectors = np.linalg.eig(u)
    scores = []
    for i, value in enumerate(values):
        p = np.abs(vectors[:, i]) ** 2
        edge_weight = float(p[:2].sum() + p[-2:].sum())
        phase_distance = abs(float(np.angle(value)))
        scores.append((edge_weight, -phase_distance, i))
    index = max(scores)[2]
    vector = vectors[:, index].astype(complex)
    vector /= np.linalg.norm(vector)
    pivot = int(np.argmax(np.abs(vector)))
    vector *= np.exp(-1j * np.angle(vector[pivot]))
    return vector


def initial_state(n_cells, n_max, mode):
    dim = n_max + 1
    n_modes = 2 * n_cells
    if mode == "edge_site":
        amplitudes = np.zeros(n_modes, dtype=complex)
        amplitudes[0] = ALPHA
    elif mode == "topological_zero_mode":
        amplitudes = ALPHA * edge_mode_vector(n_cells, CANDIDATES["chiral_topological"])
    else:
        raise ValueError(mode)
    plus = qt.tensor(*[qt.coherent(dim, value) for value in amplitudes])
    minus = qt.tensor(*[qt.coherent(dim, -value) for value in amplitudes])
    return (plus + minus).unit()


def rate_vector(n_modes, profile, disorder, seed, gamma_edge, gamma_bulk):
    base = np.full(n_modes, gamma_bulk, dtype=float)
    if profile == "edge_low_loss":
        base[[0, -1]] = gamma_edge
    elif profile != "uniform_loss":
        raise ValueError(profile)
    rng = np.random.default_rng(seed)
    return np.maximum(base + rng.uniform(-disorder, disorder) * base, 0.0)


def negativity(state, site, grid_points):
    rho = state.ptrace(site)
    if rho.isket:
        rho = qt.ket2dm(rho)
    grid = np.linspace(-4.0, 4.0, grid_points)
    w = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    return float(np.sum(np.abs(w[w < 0])) * dx * dx)


# --- Exact <=2-quanta subspace (see module docstring) -----------------------

Q_MAX = 2  # even cat on site 0 -> total quanta even, drive conserves, loss lowers


def q_max_for(n_max):
    """Minimal invariant total-quanta cutoff for the even cat initial state.

    The cat (|alpha> + |-alpha>) has even parity only, so the largest
    populated Fock number is the largest even integer <= n_max; the
    drive conserves total quanta and loss only lowers it, so Q <= that
    value is exactly invariant.  n_max=2 -> 2 (production), 3 -> 2,
    4 -> 4.

    Args:
        n_max (int): Per-site Fock truncation.

    Returns:
        int: Total-quanta subspace cutoff q_max (dimensionless).
    """
    return n_max if n_max % 2 == 0 else n_max - 1


def sub_basis(n_modes, n_max, q_max=Q_MAX):
    """Enumerate occupations with total quanta <= q_max.

    Args:
        n_modes (int): Number of modes n_modes = 2 * n_cells.
        n_max (int): Per-site Fock truncation n_max.
        q_max (int): Maximum total quanta kept.

    Returns:
        tuple[list[tuple[int, ...]], np.ndarray]: Occupation tuples (site 0
        first, matching qt.tensor ordering) and their full-space indices
        (big-endian, site 0 most significant), sorted ascending.
    """
    dim = n_max + 1
    tuples, idx = [], []
    for occ in itertools.product(range(dim), repeat=n_modes):
        if sum(occ) <= q_max:
            tuples.append(occ)
            v = 0
            for n in occ:
                v = v * dim + n
            idx.append(v)
    return tuples, np.array(idx, dtype=np.int64)


def sub_destroy(n_modes, n_max, tuples, site):
    """Annihilation operator of one site restricted to the subspace.

    Args:
        n_modes (int): Number of modes.
        n_max (int): Per-site Fock truncation.
        tuples (list[tuple[int, ...]]): Subspace occupation tuples.
        site (int): Mode index of the annihilation operator.

    Returns:
        qutip.Qobj: (dim x dim) subspace annihilation operator (dim = len(tuples)).
    """
    lookup = {t: i for i, t in enumerate(tuples)}
    rows, cols, data = [], [], []
    for j, t in enumerate(tuples):
        if t[site] > 0:
            t2 = list(t)
            amp = float(np.sqrt(t[site]))
            t2[site] -= 1
            i = lookup[tuple(t2)]
            rows.append(i)
            cols.append(j)
            data.append(amp)
    mat = sp.csr_matrix((data, (rows, cols)), shape=(len(tuples), len(tuples)))
    return qt.Qobj(mat)


def build_sub_ops(n_modes, n_max, q_max=Q_MAX):
    """Build the subspace annihilation operators for all modes.

    Args:
        n_modes (int): Number of modes.
        n_max (int): Per-site Fock truncation.
        q_max (int): Total-quanta subspace cutoff (see :func:`q_max_for`).

    Returns:
        tuple[list[qutip.Qobj], np.ndarray]: Subspace destroy operators and the
        full-space index array of the subspace (for state embedding).
    """
    tuples, idx = sub_basis(n_modes, n_max, q_max=q_max)
    return [sub_destroy(n_modes, n_max, tuples, site) for site in range(n_modes)], idx


def embed_dm(sub_state, idx, n_modes, n_max):
    """Embed a subspace ket/density matrix back into the full tensor space.

    Args:
        sub_state (qutip.Qobj): State on the subspace (ket or density matrix).
        np.ndarray: Full-space indices of the subspace (as from sub_basis).
        n_modes (int): Number of modes.
        n_max (int): Per-site Fock truncation.

    Returns:
        qutip.Qobj: State/density matrix on dims [[n_max+1]*n_modes]*2.
    """
    dm = qt.ket2dm(sub_state) if sub_state.isket else sub_state
    coo = dm.data.as_scipy().tocoo()
    n_full = (n_max + 1) ** n_modes
    mat = sp.csr_matrix((coo.data, (idx[coo.row], idx[coo.col])),
                        shape=(n_full, n_full))
    dims = [[n_max + 1] * n_modes, [n_max + 1] * n_modes]
    return qt.Qobj(mat, dims=dims)


def simulate_single(args_tuple):
    """Run one (profile, disorder, sample, sector) simulation. Worker function for Pool.

    Dynamics run exactly in the total-quanta Q<=q_max subspace (module docstring);
    the returned row schema is identical to the full-space formulation.
    """
    (profile, disorder, sample, sector_key, n_cells, n_max, n_mc,
     ntraj, wigner_grid, gamma_edge, gamma_bulk, t_stop, n_time,
     kerr, base_seed, ops_cache_key) = args_tuple

    n_modes = 2 * n_cells
    # Subspace ops rebuilt in worker (Qobj not picklable across processes)
    q_max = q_max_for(n_max)
    ops, idx = build_sub_ops(n_modes, n_max, q_max=q_max)

    candidate = CANDIDATES[sector_key]
    psi_full = initial_state(n_cells, n_max, "edge_site").full().ravel()
    psi0 = qt.Qobj(psi_full[idx])
    leakage = float(np.sqrt(max(np.sum(np.abs(psi_full) ** 2) - float(np.sum(np.abs(psi0.full().ravel()) ** 2)), 0.0)))
    assert leakage < 1e-12, f"initial state leaves the Q<={q_max} subspace: {leakage:.3e}"

    profile_list = list(DISORDER_LEVELS)
    profile_index = ("edge_low_loss", "uniform_loss").index(profile)
    disorder_index = profile_list.index(disorder)

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
    sub_state = result.states[-1]
    occ = np.asarray([float(np.real(qt.expect(op.dag() * op, sub_state))) for op in ops])

    state = embed_dm(sub_state, idx, n_modes, n_max)
    neg_edge = negativity(state, 0, wigner_grid)
    neg_center = negativity(state, n_modes // 2, wigner_grid)

    row = {
        "n_cells": n_cells, "n_modes": n_modes, "n_max": n_max,
        "sector": sector_key, **candidate, "kerr": kerr,
        "loss_profile": profile, "disorder_level": disorder,
        "sample": sample, "draw_seed": draw_seed,
        "trajectory_seed": seed, "ntraj": ntraj,
        "t_stop": t_stop, "period": PERIOD,
        "wigner_grid": wigner_grid, "rates": rates.tolist(),
        "negativity_edge": neg_edge,
        "negativity_center": neg_center,
        "edge_population": float(occ[0] + occ[-1]),
        "bulk_population": float(np.sum(occ[1:-1])),
        "occupations": occ.tolist(),
    }
    return row


def simulate_indexed(item):
    """Wrap simulate_single with the task index for checkpoint bookkeeping."""
    idx, task = item
    return idx, simulate_single(task)


def main():
    p = argparse.ArgumentParser(description="Parallel 8-mode NMC=500 chiral Floquet dynamics")
    p.add_argument("--n-cells", type=int, default=4)
    p.add_argument("--n-max", type=int, default=2)
    p.add_argument("--n-mc", type=int, default=500)
    p.add_argument("--ntraj", type=int, default=10)
    p.add_argument("--wigner-grid", type=int, default=31)
    p.add_argument("--kerr", type=float, default=0.5)
    p.add_argument("--gamma-edge", type=float, default=1e-4)
    p.add_argument("--gamma-bulk", type=float, default=0.2)
    p.add_argument("--t-stop", type=float, default=PERIOD)
    p.add_argument("--n-time", type=int, default=9)
    p.add_argument("--seed", type=int, default=20260814)
    p.add_argument("--n-workers", type=int, default=16)
    p.add_argument("--output", type=Path, default=ROOT / "data" / "chiral_production" / "chiral_floquet_4cells_NMC500.json")
    args = p.parse_args()

    out = Path(args.output)
    if out.exists():
        print(f"output already finalized ({out}); nothing to do")
        return

    n_modes = 2 * args.n_cells
    n_tasks = 2 * len(DISORDER_LEVELS) * args.n_mc * 2  # profiles × disorders × samples × sectors
    print(f"Total tasks: {n_tasks}  |  Workers: {args.n_workers}  |  n_cells: {args.n_cells}  |  NMC: {args.n_mc}")

    # Fail-closed: the Q<=q_max subspace must be exactly invariant for this design
    # (Hamiltonian block-diagonal, loss maps the subspace into itself).
    q_max = q_max_for(args.n_max)
    full_ops = build_ops(n_modes, args.n_max)
    _, sub_idx = build_sub_ops(n_modes, args.n_max, q_max=q_max)
    inside = np.zeros((args.n_max + 1) ** n_modes, dtype=bool)
    inside[sub_idx] = True
    outside = ~inside
    for cand in CANDIDATES.values():
        for h_full in (step_hamiltonian(full_ops, cand["intra_A"], 0.0, args.kerr),
                       step_hamiltonian(full_ops, 0.0, cand["inter_B"], args.kerr)):
            csr = h_full.data.as_scipy().tocsr()
            assert csr[outside][:, sub_idx].nnz == 0, "Hamiltonian leaves the Q subspace"
    for op in full_ops:
        assert op.data.as_scipy().tocsr()[outside][:, sub_idx].nnz == 0, "loss leaves the Q subspace"
    print(f"invariance OK: subspace dim {len(sub_idx)} of {(args.n_max + 1) ** n_modes}")

    # Build task list
    tasks = []
    for profile in ("edge_low_loss", "uniform_loss"):
        for disorder in DISORDER_LEVELS:
            for sample in range(args.n_mc):
                for sector_key in CANDIDATES:
                    tasks.append((
                        profile, disorder, sample, sector_key,
                        args.n_cells, args.n_max, args.n_mc,
                        args.ntraj, args.wigner_grid,
                        args.gamma_edge, args.gamma_bulk,
                        args.t_stop, args.n_time,
                        args.kerr, args.seed, None,
                    ))

    start = time.time()
    out.parent.mkdir(parents=True, exist_ok=True)
    # Row-level checkpoint: one JSON line per finished task; resume skips done idx.
    partial = Path(str(out) + ".partial.jsonl")
    done = {}
    if partial.exists():
        with open(partial) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    rec = json.loads(line)
                    done[rec["idx"]] = rec["row"]
        print(f"resume: {len(done)}/{n_tasks} rows already in {partial}")

    pending = [(i, t) for i, t in enumerate(tasks) if i not in done]
    print(f"pending: {len(pending)} rows")
    with open(partial, "a") as ck:
        with Pool(processes=args.n_workers) as pool:
            for j, (idx, row) in enumerate(
                    pool.imap_unordered(simulate_indexed, pending)):
                done[idx] = row
                ck.write(json.dumps({"idx": idx, "row": row}) + "\n")
                ck.flush()
                if (j + 1) % 200 == 0:
                    elapsed = time.time() - start
                    rate = (j + 1) / max(elapsed, 1e-9)
                    eta = (len(pending) - j - 1) / rate / 60
                    print(f"  [{j+1}/{len(pending)}] {rate:.1f} tasks/s  ETA {eta:.0f} min")

    assert len(done) == n_tasks, f"checkpoint incomplete: {len(done)}/{n_tasks}"
    rows = [done[i] for i in sorted(done)]  # canonical task order
    elapsed = time.time() - start
    print(f"Completed {len(rows)} tasks in {elapsed:.0f}s ({len(rows)/elapsed:.1f} tasks/s)")

    output = {
        "status": "COMPLETED_CHIRAL_FLOQUET_DYNAMICS_PILOT",
        "claim_scope": "finite-time dynamics in an independently audited chiral Floquet SSH model",
        "physical_sector_certified": True,
        "not_supported": ["topological protection theorem", "archived-trace validation", "experimental validation", "metrological gain"],
        "audit_source": "data/chiral_floquet_phase_audit.json",
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "design": {
            "n_cells": args.n_cells, "n_modes": n_modes, "n_max": args.n_max,
            "initial_state_mode": "edge_site",
            "kerr": args.kerr,
            "n_mc": args.n_mc, "ntraj": args.ntraj, "wigner_grid": args.wigner_grid,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": args.t_stop, "n_time_points": args.n_time,
            "period": PERIOD, "paired_disorder": True,
            "independent_unit": "disorder realization; trajectories nested",
            "disorder_levels": list(DISORDER_LEVELS),
            "n_workers": args.n_workers,
            "state_reduction": (
                f"exact total-quanta Q<={q_max} subspace (dim {len(sub_idx)} "
                f"of {(args.n_max + 1) ** n_modes}); invariance checked at run start; "
                "record-level validated against full-space rows 2026-10-02 (diffs <= 3e-5)"
            ),
        },
        "records": rows,
        "elapsed_seconds": elapsed,
    }
    out.write_text(json.dumps(output, indent=2) + "\n")
    partial.unlink()  # checkpoint consumed by the finalized output
    print(f"saved {out} with {len(rows)} records")


if __name__ == "__main__":
    main()
