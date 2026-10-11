#!/usr/bin/env python3
"""Dissipative dynamics for the independently audited chiral Floquet SSH model.

This model is separate from the archived traces and from the legacy FL_QOM
massive two-step model.  Its bulk invariants and finite-chain edge diagnostics
are recorded in data/chiral_floquet_phase_audit.json.
"""
from __future__ import annotations

import argparse
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
    # Deterministic global phase for provenance and reproducibility.
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
    """Sample per-site energy-loss rates for one disorder draw.

    Args:
        n_modes (int): Number of lattice sites (dimensionless).
        profile (str): Loss profile. ``edge_low_loss`` damps all sites at
            ``gamma_bulk`` except the two edge sites (``gamma_edge``);
            ``uniform_loss`` damps every site at ``gamma_bulk``;
            ``chiral_even_loss`` damps only the even (A) sublattice
            (site 0 at ``gamma_edge``, other even sites at ``gamma_bulk``)
            so that every collapse operator $L_k$ commutes with the chiral
            operator $S$ ($S L_k S = L_k$); ``chiral_odd_loss`` damps only
            the odd (B) sublattice at ``gamma_bulk``;
            ``chiral_even_loss_matched`` is ``chiral_even_loss`` rescaled by
            $n_{\rm modes}/n_{\rm damped} = 2$, so the total loss-rate budget
            and the edge/bulk rate ratio match ``edge_low_loss`` exactly.
        disorder (float): Relative multiplicative noise amplitude
            (dimensionless, in units of the site rate).
        seed (int): Deterministic seed of the disorder draw (dimensionless).
        gamma_edge (float): Edge-site loss rate (units of drive frequency).
        gamma_bulk (float): Bulk-site loss rate (units of drive frequency).

    Returns:
        numpy.ndarray: Non-negative loss rates per site, shape
        ``(n_modes,)``, units of drive frequency.
    """
    base = np.zeros(n_modes, dtype=float)
    if profile == "edge_low_loss":
        base[:] = gamma_bulk
        base[[0, -1]] = gamma_edge
    elif profile == "uniform_loss":
        base[:] = gamma_bulk
    elif profile == "chiral_even_loss":
        base[0::2] = gamma_bulk
        base[0] = gamma_edge
    elif profile == "chiral_odd_loss":
        base[1::2] = gamma_bulk
    elif profile == "chiral_even_loss_matched":
        base[0::2] = 2.0 * gamma_bulk
        base[0] = 2.0 * gamma_edge
    else:
        raise ValueError(profile)
    rng = np.random.default_rng(seed)
    return np.maximum(base + rng.uniform(-disorder, disorder) * base, 0.0)


def negativity(state, site, grid_points):
    """Integrated absolute Wigner negativity of one site of a bipartite state.

    Args:
        state (qutip.Qobj): Density matrix (or ket) of the full $n$-mode state.
        site (int): Site index whose reduced state is evaluated (dimensionless).
        grid_points (int): Side length of the square Wigner grid (dimensionless).

    Returns:
        float: $\int_{w<0} |w(x,p)|\,dx\,dp$ of the reduced Wigner function,
        dimensionless (grid spans $[-4, 4]$ in $\hbar$-scaled units).
    """
    rho = state.ptrace(site)
    if rho.isket:
        rho = qt.ket2dm(rho)
    grid = np.linspace(-4.0, 4.0, grid_points)
    w = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    return float(np.sum(np.abs(w[w < 0])) * dx * dx)


def simulate(h, psi0, c_ops, ops, ntraj, seed, t_stop, n_time, grid_points,
             series=False):
    """Propagate the stochastic master equation and evaluate the observables.

    Args:
        h (qutip.QobjEvo): Time-dependent Floquet Hamiltonian (units of drive frequency).
        psi0 (qutip.Qobj): Initial even-cat state (ket).
        c_ops (list): Collapse operators $\sqrt{\gamma_j} a_j$ (units of drive frequency).
        ops (list): Single-site annihilation operators of the truncated Fock space.
        ntraj (int): Number of Monte Carlo trajectories (dimensionless).
        seed (int): Trajectory seed (dimensionless).
        t_stop (float): Final evolution time (units of inverse drive frequency, i.e. periods/$2\pi$).
        n_time (int): Number of stored evaluation times in $[0, t_{\rm stop}]$ (dimensionless).
        grid_points (int): Side length of the Wigner grid (dimensionless).
        series (bool): If True, also return the observables at every stored
            evaluation time (micromotion scan).

    Returns:
        dict: Terminal-time observables (Wigner-void negativity of the edge
        and centre mode, edge/bulk populations, occupations) and, when
        ``series`` is True, a ``time_series`` block with the same
        observables evaluated at each stored time.
    """
    times = np.linspace(0.0, t_stop, n_time)
    result = qt.mcsolve(
        h, psi0, times, c_ops,
        ntraj=ntraj, seeds=seed,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    states = result.states
    state = states[-1]
    occ = np.asarray([float(np.real(qt.expect(op.dag() * op, state))) for op in ops])
    out = {
        "negativity_edge": negativity(state, 0, grid_points),
        "negativity_center": negativity(state, len(ops) // 2, grid_points),
        "edge_population": float(occ[0] + occ[-1]),
        "bulk_population": float(np.sum(occ[1:-1])),
        "occupations": occ.tolist(),
    }
    if series:
        center = len(ops) // 2
        edge_pop, bulk_pop, neg_edge, neg_center = [], [], [], []
        for st in states:
            st_occ = np.asarray([float(np.real(qt.expect(op.dag() * op, st))) for op in ops])
            edge_pop.append(float(st_occ[0] + st_occ[-1]))
            bulk_pop.append(float(np.sum(st_occ[1:-1])))
            neg_edge.append(negativity(st, 0, grid_points))
            neg_center.append(negativity(st, center, grid_points))
        out["time_series"] = {
            "times": [float(t) for t in times],
            "fraction_of_period": [float(t / PERIOD) for t in times],
            "negativity_edge": neg_edge,
            "negativity_center": neg_center,
            "edge_population": edge_pop,
            "bulk_population": bulk_pop,
        }
    return out


def run(args):
    n_modes = 2 * args.n_cells
    ops = build_ops(n_modes, args.n_max)
    psi0 = initial_state(args.n_cells, args.n_max, args.initial_state_mode)
    rows = []
    start = time.time()
    for profile_index, profile in enumerate(("edge_low_loss", "uniform_loss")):
        for disorder_index, disorder in enumerate(DISORDER_LEVELS):
            for sample in range(args.n_mc):
                draw_seed = args.seed + profile_index * 1_000_000 + disorder_index * 10_000 + sample
                rates = rate_vector(n_modes, profile, disorder, draw_seed,
                                   args.gamma_edge, args.gamma_bulk)
                c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
                for sector_index, (sector, candidate) in enumerate(CANDIDATES.items()):
                    seed = draw_seed + 100 + sector_index
                    row = simulate(
                        driven_hamiltonian(ops, candidate, args.kerr), psi0, c_ops, ops,
                        args.ntraj, seed, args.t_stop, args.n_time, args.wigner_grid,
                    )
                    rows.append({
                        "n_cells": args.n_cells, "n_modes": n_modes, "n_max": args.n_max,
                        "sector": sector, **candidate, "kerr": args.kerr,
                        "loss_profile": profile, "disorder_level": disorder,
                        "sample": sample, "draw_seed": draw_seed,
                        "trajectory_seed": seed, "ntraj": args.ntraj,
                        "t_stop": args.t_stop, "period": PERIOD,
                        "wigner_grid": args.wigner_grid, "rates": rates.tolist(), **row,
                    })
            print(f"completed profile={profile} disorder={disorder}")
    output = {
        "status": "COMPLETED_CHIRAL_FLOQUET_DYNAMICS_PILOT",
        "claim_scope": "finite-time dynamics in an independently audited chiral Floquet SSH model",
        "physical_sector_certified": True,
        "not_supported": ["topological protection theorem", "archived-trace validation", "experimental validation", "metrological gain"],
        "audit_source": "data/chiral_floquet_phase_audit.json",
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "design": {
            "n_cells": args.n_cells, "n_modes": n_modes, "n_max": args.n_max,
            "initial_state_mode": args.initial_state_mode,
            "kerr": args.kerr,
            "n_mc": args.n_mc, "ntraj": args.ntraj, "wigner_grid": args.wigner_grid,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": args.t_stop, "n_time_points": args.n_time,
            "period": PERIOD, "paired_disorder": True,
            "independent_unit": "disorder realization; trajectories nested",
            "disorder_levels": DISORDER_LEVELS,
        },
        "records": rows,
        "elapsed_seconds": time.time() - start,
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(output, indent=2) + "\n")
    print(f"saved {out} with {len(rows)} records")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-cells", type=int, default=2)
    p.add_argument("--n-max", type=int, default=2)
    p.add_argument("--n-mc", type=int, default=2)
    p.add_argument("--ntraj", type=int, default=5)
    p.add_argument("--wigner-grid", type=int, default=31)
    p.add_argument("--initial-state-mode", choices=("edge_site", "topological_zero_mode"), default="edge_site")
    p.add_argument("--kerr", type=float, default=0.5)
    p.add_argument("--gamma-edge", type=float, default=1e-4)
    p.add_argument("--gamma-bulk", type=float, default=0.2)
    p.add_argument("--t-stop", type=float, default=PERIOD)
    p.add_argument("--n-time", type=int, default=9)
    p.add_argument("--seed", type=int, default=20260814)
    p.add_argument("--output", type=Path, default=ROOT / "data" / "chiral_floquet_pilot.json")
    args = p.parse_args()
    if min(args.n_cells, args.n_max, args.n_mc, args.ntraj) < 1 or args.kerr < 0 or args.t_stop <= 0 or args.n_time < 2:
        p.error("invalid simulation arguments")
    run(args)


if __name__ == "__main__":
    main()
