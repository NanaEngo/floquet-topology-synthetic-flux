#!/usr/bin/env python3
"""Pre-specified edge-mode mechanism pilot for the chiral Floquet SSH model.

This pilot does not claim a topological protection theorem.  It tests whether the
same certified Floquet backbone produces a measurable edge-mode survival signal
under matched dissipation, and separates the passive limit (where an effective
loss bound can be checked) from the Kerr extension.
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
ALPHA = 1.0
CANDIDATES = {
    "chiral_topological": {"intra_A": 0.1, "inter_B": 1.2, "nu_0": -1, "nu_pi": 0},
    "chiral_trivial": {"intra_A": 1.2, "inter_B": 0.1, "nu_0": 0, "nu_pi": 0},
}
DISORDER_LEVELS = (0.0, 0.5, 1.0, 1.5)


def finite_chain_matrix(n_cells: int, intra: float, inter: float) -> np.ndarray:
    n = 2 * n_cells
    h = np.zeros((n, n), dtype=complex)
    for cell in range(n_cells):
        a, b = 2 * cell, 2 * cell + 1
        h[a, b] = h[b, a] = intra
        if cell < n_cells - 1:
            an, bn = 2 * (cell + 1), 2 * (cell + 1) + 1
            h[a, bn] = h[bn, a] = inter
    return h


def edge_mode_vector(n_cells: int, candidate: dict) -> np.ndarray:
    """Return the most edge-localized Floquet eigenvector for the certified sector."""
    ha = finite_chain_matrix(n_cells, candidate["intra_A"], 0.0)
    hb = finite_chain_matrix(n_cells, 0.0, candidate["inter_B"])
    u = expm(-1j * hb * HALF_PERIOD) @ expm(-1j * ha * HALF_PERIOD)
    values, vectors = np.linalg.eig(u)
    scores = []
    for i, value in enumerate(values):
        p = np.abs(vectors[:, i]) ** 2
        edge_weight = float(p[:2].sum() + p[-2:].sum())
        scores.append((edge_weight, -abs(float(np.angle(value))), i))
    vector = vectors[:, max(scores)[2]].astype(complex)
    vector /= np.linalg.norm(vector)
    pivot = int(np.argmax(np.abs(vector)))
    vector *= np.exp(-1j * np.angle(vector[pivot]))
    return vector


def mode_number_operator(vector: np.ndarray, ops, convention: str = "standard"):
    """Return b^†b for b=sum_j conjugate(vector[j]) a_j.

    ``standard`` is the mathematically correct coefficient convention.  The
    ``legacy_fixed`` option is retained only to reproduce the earlier fixed-frame
    pilot outputs; it must not be used for the covariant observable.
    """
    if convention == "standard":
        coefficient = lambda i, j: vector[i] * np.conjugate(vector[j])
    elif convention == "legacy_fixed":
        coefficient = lambda i, j: np.conjugate(vector[i]) * vector[j]
    else:
        raise ValueError(f"unknown edge operator convention: {convention}")
    operator = 0 * ops[0]
    for i in range(len(ops)):
        for j in range(len(ops)):
            operator += coefficient(i, j) * ops[i].dag() * ops[j]
    return operator


def micromotion_vector(t: float, n_cells: int, candidate: dict, vector: np.ndarray) -> np.ndarray:
    """Transport a t=0 Floquet eigenvector through the coherent micromotion."""
    tau = float(t) % PERIOD
    ha = finite_chain_matrix(n_cells, candidate["intra_A"], 0.0)
    hb = finite_chain_matrix(n_cells, 0.0, candidate["inter_B"])
    if tau < HALF_PERIOD:
        propagator = expm(-1j * ha * tau)
    else:
        propagator = expm(-1j * hb * (tau - HALF_PERIOD)) @ expm(-1j * ha * HALF_PERIOD)
    transported = propagator @ vector
    return transported / np.linalg.norm(transported)


def build_ops(n_modes: int, n_max: int):
    a = qt.destroy(n_max + 1)
    eye = qt.qeye(n_max + 1)
    return [qt.tensor(*[a if i == j else eye for i in range(n_modes)])
            for j in range(n_modes)]


def second_quantize(matrix: np.ndarray, ops):
    h = 0 * ops[0]
    for i in range(matrix.shape[0]):
        for j in range(i + 1, matrix.shape[1]):
            value = matrix[i, j]
            if abs(value) > 0:
                h += value * ops[i].dag() * ops[j]
                h += np.conjugate(value) * ops[j].dag() * ops[i]
    return h


def step_hamiltonian(ops, intra: float, inter: float, kerr: float):
    h = second_quantize(finite_chain_matrix(len(ops) // 2, intra, inter), ops)
    # The Kerr term is deliberately local and is recorded as a boundary
    # perturbation; it is not part of the chiral single-particle invariant.
    h += kerr * ops[0].dag() * ops[0].dag() * ops[0] * ops[0]
    return h


def coeff_a(t, **kwargs):
    return 1.0 if (float(t) % PERIOD) < HALF_PERIOD else 0.0


def coeff_b(t, **kwargs):
    return 1.0 - coeff_a(t, **kwargs)


def driven_hamiltonian(ops, candidate: dict, kerr: float):
    ha = step_hamiltonian(ops, candidate["intra_A"], 0.0, kerr)
    hb = step_hamiltonian(ops, 0.0, candidate["inter_B"], kerr)
    return qt.QobjEvo([[ha, coeff_a], [hb, coeff_b]])


def cat_state(n_cells: int, n_max: int, preparation: str, edge_vector: np.ndarray):
    n_modes = 2 * n_cells
    if preparation == "common_edge_site":
        amplitudes = np.zeros(n_modes, dtype=complex)
        amplitudes[0] = ALPHA
    elif preparation == "topological_edge_mode":
        amplitudes = ALPHA * edge_vector
    else:
        raise ValueError(preparation)
    dim = n_max + 1
    plus = qt.tensor(*[qt.coherent(dim, value) for value in amplitudes])
    minus = qt.tensor(*[qt.coherent(dim, -value) for value in amplitudes])
    return (plus + minus).unit()


def rate_vector(n_modes: int, profile: str, disorder: float, seed: int,
                gamma_edge: float, gamma_bulk: float) -> np.ndarray:
    base = np.full(n_modes, gamma_bulk, dtype=float)
    if profile == "edge_low_loss":
        base[[0, -1]] = gamma_edge
    elif profile != "uniform_loss":
        raise ValueError(profile)
    rng = np.random.default_rng(seed)
    return np.maximum(base + rng.uniform(-disorder, disorder, size=n_modes) * base, 0.0)


def expect_real(operator, state) -> float:
    return float(np.real(qt.expect(operator, state)))


def local_wigner_negativity(state, site: int, grid_points: int) -> float:
    rho = state.ptrace(site)
    if rho.isket:
        rho = qt.ket2dm(rho)
    grid = np.linspace(-4.0, 4.0, grid_points)
    w = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    return float(np.sum(np.abs(w[w < 0])) * dx * dx)


def simulate(h, psi0, c_ops, ops, edge_operator_series, rates, ntraj, seed,
             t_stop, n_time, grid_points, passive_prediction_series, passive_status):
    result = qt.mcsolve(
        h, psi0, np.linspace(0.0, t_stop, n_time), c_ops,
        ntraj=ntraj, seeds=seed,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    state = result.states[-1]
    times = np.linspace(0.0, t_stop, n_time)
    edge_population_series = [expect_real(operator, state_i)
                              for operator, state_i in zip(edge_operator_series, result.states)]
    passive_error_series = (None if passive_prediction_series is None else
                            [float(a - b) for a, b in zip(edge_population_series, passive_prediction_series)])
    occupations = np.asarray([expect_real(op.dag() * op, state) for op in ops])
    edge_population = edge_population_series[-1]
    total_population = float(np.sum(occupations))
    edge_loss = float(rates[0] * occupations[0] + rates[-1] * occupations[-1])
    bulk_loss = float(np.sum(rates[1:-1] * occupations[1:-1]))
    return {
        "edge_mode_population": edge_population,
        "total_population": total_population,
        "bulk_mode_population_proxy": max(total_population - edge_population, 0.0),
        "edge_loss_current": edge_loss,
        "bulk_loss_current": bulk_loss,
        "negativity_site_0": local_wigner_negativity(state, 0, grid_points),
        "negativity_site_center": local_wigner_negativity(state, len(ops) // 2, grid_points),
        "occupations": occupations.tolist(),
        "time_points": times.tolist(),
        "edge_mode_population_timeseries": edge_population_series,
        "passive_prediction_timeseries": (None if passive_prediction_series is None else list(passive_prediction_series)),
        "passive_prediction_edge_population": (None if passive_prediction_series is None else passive_prediction_series[-1]),
        "passive_prediction_error": (None if passive_error_series is None else passive_error_series[-1]),
        "max_abs_passive_residual": (None if passive_error_series is None else float(max(abs(x) for x in passive_error_series))),
        "passive_reference_status": passive_status,
    }


def run(args):
    n_modes = 2 * args.n_cells
    ops = build_ops(n_modes, args.n_max)
    topological_candidate = CANDIDATES["chiral_topological"]
    edge_vector = edge_mode_vector(args.n_cells, topological_candidate)
    psi0 = cat_state(args.n_cells, args.n_max, args.initial_state_mode, edge_vector)
    rows = []
    start = time.time()
    times = np.linspace(0.0, args.t_stop, args.n_time)
    if args.edge_observable_frame == "micromotion":
        edge_vectors = [micromotion_vector(t, args.n_cells, topological_candidate, edge_vector)
                        for t in times]
        operator_convention = "standard"
    else:
        edge_vectors = [edge_vector] * len(times)
        operator_convention = args.edge_operator_convention
    edge_operator_series = [mode_number_operator(v, ops, operator_convention)
                            for v in edge_vectors]
    initial_edge_population = expect_real(edge_operator_series[0], psi0)
    for profile_index, profile in enumerate(args.profiles):
        for disorder_index, disorder in enumerate(args.disorder_levels):
            for sample in range(args.n_mc):
                draw_seed = args.seed + profile_index * 1_000_000 + disorder_index * 10_000 + sample
                rates = rate_vector(n_modes, profile, disorder, draw_seed,
                                   args.gamma_edge, args.gamma_bulk)
                c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
                gamma_eff_series = np.asarray(
                    [np.sum(rates * np.abs(v) ** 2) for v in edge_vectors], dtype=float
                )
                integrated_exposure = np.zeros_like(times)
                integrated_exposure[1:] = np.cumsum(
                    0.5 * (gamma_eff_series[:-1] + gamma_eff_series[1:]) * np.diff(times)
                )
                gamma_eff = float(gamma_eff_series[0])
                prediction_series = initial_edge_population * np.exp(-integrated_exposure)
                for sector_index, (sector, candidate) in enumerate(CANDIDATES.items()):
                    trajectory_seed = draw_seed + 100 + sector_index
                    passive_applicable = sector == "chiral_topological" and args.kerr == 0.0
                    record = simulate(
                        driven_hamiltonian(ops, candidate, args.kerr), psi0, c_ops,
                        ops, edge_operator_series, rates, args.ntraj, trajectory_seed,
                        args.t_stop, args.n_time, args.wigner_grid,
                        prediction_series if passive_applicable else None,
                        "topological_passive_reference" if passive_applicable else "comparison_only",
                    )
                    rows.append({
                        "n_cells": args.n_cells, "n_modes": n_modes, "n_max": args.n_max,
                        "sector": sector, **candidate, "kerr": args.kerr,
                        "initial_state_mode": args.initial_state_mode,
                        "loss_profile": profile, "disorder_level": disorder,
                        "sample": sample, "draw_seed": draw_seed,
                        "trajectory_seed": trajectory_seed, "ntraj": args.ntraj,
                        "t_stop": args.t_stop, "period": PERIOD,                        "wigner_grid": args.wigner_grid, "rates": rates.tolist(),
                        "edge_vector": [[float(value.real), float(value.imag)] for value in edge_vector],
                        "edge_observable_frame": args.edge_observable_frame,
                        "edge_operator_convention": operator_convention,
                        "gamma_eff": gamma_eff,
                        "gamma_eff_timeseries": gamma_eff_series.tolist(),

                        "initial_edge_population": initial_edge_population, **record,
                    })
            print(f"completed preparation={args.initial_state_mode} kerr={args.kerr} profile={profile} disorder={disorder}")
    output = {
        "status": "COMPLETED_EDGE_MODE_MECHANISM_PILOT",
        "scientific_status": "PROVISIONAL_EDGE_MODE_OBSERVABLE_PILOT",
        "claim_scope": "finite-time matched dynamics; no topological-protection or metrological claim",
        "model": "same two-step chiral Floquet SSH backbone as the invariant audit, with local Kerr boundary perturbation and site Lindblad losses",
        "audit_source": "data/chiral_floquet_phase_audit.json",
        "parameter_status": "NORMALIZED_UNCALIBRATED_NOT_OPTIMIZED",
        "edge_observable_frame": args.edge_observable_frame,
        "edge_operator_convention": operator_convention,
        "reference_frequency": None,
        "physical_units": None,
        "optimality_status": "NOT_OPTIMIZED; sectors selected for invariant/gap diagnostics, not an objective optimum",
        "design": {
            "n_cells": args.n_cells, "n_modes": n_modes, "n_max": args.n_max,
            "n_mc": args.n_mc, "ntraj": args.ntraj, "wigner_grid": args.wigner_grid,
            "initial_state_mode": args.initial_state_mode, "kerr": args.kerr,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": args.t_stop, "n_time_points": args.n_time,
            "period": PERIOD, "paired_disorder": True,
            "profiles": list(args.profiles),
            "disorder_levels": list(args.disorder_levels),
            "independent_unit": "paired disorder realization; trajectories nested",
            "passive_bound": "n_edge(t) = n_edge(0) exp(-gamma_eff t) is recorded only for the topological sector at kerr=0; it is a reference diagnostic, exact only for an isolated passive edge mode",
        },
        "edge_mode_definition": {
            "operator": "b_edge = sum_j conjugate(u_j) a_j",
            "vector_source": "most edge-localized eigenvector of the certified topological one-period single-particle unitary",
            "local_wigner_status": "site-resolved Wigner negativity is recorded; collective-mode Wigner negativity requires a validated bosonic mode transform and is not silently substituted",
        },
        "records": rows,
        "elapsed_seconds": time.time() - start,
        "not_supported": [
            "global basin structure", "topological protection theorem", "Liouvillian invariant",
            "experimental validation", "metrological gain", "collective-mode Wigner negativity until transform validation",
        ],
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
    p.add_argument("--ntraj", type=int, default=2)
    p.add_argument("--wigner-grid", type=int, default=31)
    p.add_argument("--initial-state-mode", choices=("common_edge_site", "topological_edge_mode"), default="common_edge_site")
    p.add_argument("--profiles", nargs="+", choices=("edge_low_loss", "uniform_loss"), default=["edge_low_loss", "uniform_loss"])
    p.add_argument("--disorder-levels", nargs="+", type=float, default=list(DISORDER_LEVELS))
    p.add_argument("--kerr", type=float, default=0.0)
    p.add_argument("--gamma-edge", type=float, default=1e-4)
    p.add_argument("--gamma-bulk", type=float, default=0.2)
    p.add_argument("--t-stop", type=float, default=8.0 * PERIOD)
    p.add_argument("--n-time", type=int, default=17)
    p.add_argument("--seed", type=int, default=20260817)
    p.add_argument("--edge-observable-frame", choices=("fixed", "micromotion"), default="fixed")
    p.add_argument("--edge-operator-convention", choices=("standard", "legacy_fixed"), default="legacy_fixed")
    p.add_argument("--output", type=Path, default=ROOT / "data" / "edge_mode_pilot.json")
    args = p.parse_args()
    if min(args.n_cells, args.n_max, args.n_mc, args.ntraj) < 1 or args.kerr < 0 or args.t_stop <= 0 or args.n_time < 2:
        p.error("invalid simulation arguments")
    run(args)


if __name__ == "__main__":
    main()
