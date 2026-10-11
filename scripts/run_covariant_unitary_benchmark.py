#!/usr/bin/env python3
"""Deterministic unitary benchmark for the covariant Floquet edge observable.

This test is deliberately separate from the open-system Monte Carlo pilot.  It asks
whether the time-dependent edge vector u_e(t)=U(t)u_e(0) and the operator
b_e(t)=sum_j conjugate(u_e,j(t)) a_j are mathematically covariant under the same
single-particle drive used to define the Floquet edge mode.

A PASS here validates the observable/frame correspondence in the unitary model.  It
does not validate the lossy passive reference, the Kerr dynamics, a collective-mode
Wigner transform, or topological protection.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm, null_space

ROOT = Path(__file__).resolve().parents[1]
PERIOD = 2.0 * np.pi
HALF_PERIOD = np.pi
CANDIDATE = {"intra_A": 0.1, "inter_B": 1.2}


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


def step_unitary(n_cells: int, intra: float, inter: float, duration: float) -> np.ndarray:
    return expm(-1j * finite_chain_matrix(n_cells, intra, inter) * duration)


def period_unitary(n_cells: int) -> np.ndarray:
    ua_half = step_unitary(n_cells, CANDIDATE["intra_A"], 0.0, HALF_PERIOD)
    ub_half = step_unitary(n_cells, 0.0, CANDIDATE["inter_B"], HALF_PERIOD)
    return ub_half @ ua_half


def coherent_micromotion(t: float, n_cells: int) -> np.ndarray:
    # Include the right-hand endpoint explicitly: U(T)=U_F, not U(0).
    if float(t) > 0.0 and np.isclose(float(t) % PERIOD, 0.0, atol=1e-12):
        return period_unitary(n_cells)
    tau = float(t) % PERIOD
    ua_half = step_unitary(n_cells, CANDIDATE["intra_A"], 0.0, HALF_PERIOD)
    if tau < HALF_PERIOD:
        return step_unitary(n_cells, CANDIDATE["intra_A"], 0.0, tau)
    return step_unitary(n_cells, 0.0, CANDIDATE["inter_B"], tau - HALF_PERIOD) @ ua_half


def edge_mode_vector(n_cells: int) -> np.ndarray:
    uf = period_unitary(n_cells)
    values, vectors = np.linalg.eig(uf)
    scores = []
    for i in range(len(values)):
        p = np.abs(vectors[:, i]) ** 2
        edge_weight = float(p[:2].sum() + p[-2:].sum())
        scores.append((edge_weight, -abs(float(np.angle(values[i]))), i))
    vector = vectors[:, max(scores)[2]].astype(complex)
    vector /= np.linalg.norm(vector)
    pivot = int(np.argmax(np.abs(vector)))
    vector *= np.exp(-1j * np.angle(vector[pivot]))
    return vector


def evaluate_size(n_cells: int, n_time: int) -> dict:
    n = 2 * n_cells
    u0 = edge_mode_vector(n_cells)
    bulk_basis = null_space(u0.conj()[None, :])
    uf = coherent_micromotion(PERIOD, n_cells)
    floquet_phase = np.angle(np.vdot(u0, uf @ u0))
    floquet_residual = float(np.linalg.norm(uf @ u0 - np.exp(1j * floquet_phase) * u0))

    times = np.linspace(0.0, PERIOD, n_time)
    norm_errors = []
    unitarity_errors = []
    covariant_populations = []
    fixed_frame_populations = []
    bulk_leakage = []
    for t in times:
        ut = coherent_micromotion(float(t), n_cells)
        ue_t = ut @ u0
        bulk_t = ut @ bulk_basis
        norm_errors.append(abs(np.linalg.norm(ue_t) - 1.0))
        unitarity_errors.append(float(np.linalg.norm(ut.conj().T @ ut - np.eye(n))))
        covariant_populations.append(float(abs(np.vdot(ue_t, ue_t)) ** 2))
        fixed_frame_populations.append(float(abs(np.vdot(u0, ue_t)) ** 2))
        bulk_leakage.append(float(np.max(np.abs(np.conjugate(ue_t) @ bulk_t) ** 2)))

    return {
        "n_cells": n_cells,
        "n_modes": n,
        "n_time": n_time,
        "candidate": CANDIDATE,
        "edge_weight_initial": float(np.sum(np.abs(u0[[0, 1, -2, -1]]) ** 2)),
        "floquet_phase": float(floquet_phase),
        "floquet_eigenvector_residual": floquet_residual,
        "max_edge_vector_norm_error": float(max(norm_errors)),
        "max_unitarity_residual": float(max(unitarity_errors)),
        "covariant_population_min": float(min(covariant_populations)),
        "covariant_population_max": float(max(covariant_populations)),
        "max_covariant_population_error": float(max(abs(x - 1.0) for x in covariant_populations)),
        "fixed_frame_population_min": float(min(fixed_frame_populations)),
        "fixed_frame_population_max": float(max(fixed_frame_populations)),
        "max_covariant_edge_bulk_leakage": float(max(bulk_leakage)),
        "times": times.tolist(),
        "covariant_population": covariant_populations,
        "fixed_frame_population": fixed_frame_populations,
        "edge_bulk_leakage": bulk_leakage,
        "pass": bool(
            floquet_residual < 1e-10
            and max(norm_errors) < 1e-12
            and max(unitarity_errors) < 1e-11
            and max(abs(x - 1.0) for x in covariant_populations) < 1e-11
            and max(bulk_leakage) < 1e-11
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", nargs="+", type=int, default=[2, 3, 4, 6])
    parser.add_argument("--n-time", type=int, default=257)
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "edge_mode_unitary_benchmark.json")
    args = parser.parse_args()
    if any(x < 1 for x in args.sizes) or args.n_time < 3:
        parser.error("sizes must be positive and n-time must be at least 3")
    results = [evaluate_size(n, args.n_time) for n in args.sizes]
    output = {
        "status": "COMPLETED_DETERMINISTIC_UNITARY_COVARIANCE_BENCHMARK",
        "scientific_scope": "single-particle coherent micromotion and edge-bulk orthogonality only",
        "model": "two-step chiral Floquet SSH backbone; no loss, no Kerr, no disorder",
        "observable_definition": "u_e(t)=U_coh(t)u_e(0); b_e(t)=sum_j conjugate(u_e,j(t)) a_j",
        "thresholds": {
            "floquet_eigenvector_residual": 1e-10,
            "unitarity_residual": 1e-11,
            "covariant_population_error": 1e-11,
            "edge_bulk_leakage": 1e-11,
        },
        "results": results,
        "overall_pass": bool(all(x["pass"] for x in results)),
        "claim_boundary": [
            "validates the covariant observable in the coherent unitary model",
            "does not validate the open-system passive reference",
            "does not validate Kerr dynamics or collective-mode Wigner negativity",
            "does not establish topological protection",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"overall_pass": output["overall_pass"], "sizes": [{"n_cells": x["n_cells"], "pass": x["pass"], "max_covariant_population_error": x["max_covariant_population_error"], "max_covariant_edge_bulk_leakage": x["max_covariant_edge_bulk_leakage"], "fixed_frame_population_min": x["fixed_frame_population_min"]} for x in results]}, indent=2))


if __name__ == "__main__":
    main()
