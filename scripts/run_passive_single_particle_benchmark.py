#!/usr/bin/env python3
"""Deterministic passive single-particle benchmark.

The test removes Kerr nonlinearity, Fock truncation and trajectory sampling.  It
propagates one excitation under the same driven single-particle Hamiltonian with
site losses and compares the covariant edge population with the integrated exposure
reference.  Uniform loss is an analytic control; edge-selective loss tests whether
loss and micromotion generate edge--bulk transfer.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
PERIOD = 2.0 * np.pi
HALF_PERIOD = np.pi
GAMMA_EDGE = 1e-4
GAMMA_BULK = 0.2
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


def period_unitary(n_cells: int) -> np.ndarray:
    ha = finite_chain_matrix(n_cells, CANDIDATE["intra_A"], 0.0)
    hb = finite_chain_matrix(n_cells, 0.0, CANDIDATE["inter_B"])
    return expm(-1j * hb * HALF_PERIOD) @ expm(-1j * ha * HALF_PERIOD)


def edge_mode_vector(n_cells: int) -> np.ndarray:
    values, vectors = np.linalg.eig(period_unitary(n_cells))
    scores = []
    for i in range(len(values)):
        p = np.abs(vectors[:, i]) ** 2
        edge_weight = float(p[:2].sum() + p[-2:].sum())
        scores.append((edge_weight, -abs(float(np.angle(values[i]))), i))
    u = vectors[:, max(scores)[2]].astype(complex)
    u /= np.linalg.norm(u)
    pivot = int(np.argmax(np.abs(u)))
    return u * np.exp(-1j * np.angle(u[pivot]))


def coherent_micromotion(t: float, n_cells: int) -> np.ndarray:
    if float(t) > 0.0 and np.isclose(float(t) % PERIOD, 0.0, atol=1e-12):
        return period_unitary(n_cells)
    tau = float(t) % PERIOD
    ha = finite_chain_matrix(n_cells, CANDIDATE["intra_A"], 0.0)
    hb = finite_chain_matrix(n_cells, 0.0, CANDIDATE["inter_B"])
    if tau < HALF_PERIOD:
        return expm(-1j * ha * tau)
    return expm(-1j * hb * (tau - HALF_PERIOD)) @ expm(-1j * ha * HALF_PERIOD)


def rate_vector(n_modes: int, profile: str) -> np.ndarray:
    rates = np.full(n_modes, GAMMA_BULK, dtype=float)
    if profile == "edge_low_loss":
        rates[[0, -1]] = GAMMA_EDGE
    elif profile != "uniform_loss":
        raise ValueError(profile)
    return rates


def effective_propagator(t: float, n_cells: int, rates: np.ndarray) -> np.ndarray:
    ha = finite_chain_matrix(n_cells, CANDIDATE["intra_A"], 0.0) - 0.5j * np.diag(rates)
    hb = finite_chain_matrix(n_cells, 0.0, CANDIDATE["inter_B"]) - 0.5j * np.diag(rates)
    if float(t) > 0.0 and np.isclose(float(t) % PERIOD, 0.0, atol=1e-12):
        return expm(-1j * hb * HALF_PERIOD) @ expm(-1j * ha * HALF_PERIOD)
    tau = float(t) % PERIOD
    if tau < HALF_PERIOD:
        return expm(-1j * ha * tau)
    return expm(-1j * hb * (tau - HALF_PERIOD)) @ expm(-1j * ha * HALF_PERIOD)


def evaluate(n_cells: int, profile: str, n_time: int) -> dict:
    n = 2 * n_cells
    rates = rate_vector(n, profile)
    u0 = edge_mode_vector(n_cells)
    times = np.linspace(0.0, PERIOD, n_time)
    gamma_eff = []
    edge_population = []
    no_jump_norm = []
    bulk_population = []
    unitarity_norm = []
    for t in times:
        uc = coherent_micromotion(float(t), n_cells)
        ut = uc @ u0
        gamma_eff.append(float(np.sum(rates * np.abs(ut) ** 2)))
        K = effective_propagator(float(t), n_cells, rates)
        psi = K @ u0
        rho = np.outer(psi, np.conjugate(psi))
        edge_population.append(float(np.real(np.vdot(ut, rho @ ut))))
        no_jump_norm.append(float(np.real(np.trace(rho))))
        projector_bulk = np.eye(n) - np.outer(ut, np.conjugate(ut))
        bulk_population.append(float(np.real(np.trace(projector_bulk @ rho))))
        unitarity_norm.append(float(np.linalg.norm(uc.conj().T @ uc - np.eye(n))))
    gamma_eff = np.asarray(gamma_eff)
    times_array = np.asarray(times)
    exposure = np.zeros_like(times_array)
    exposure[1:] = np.cumsum(0.5 * (gamma_eff[:-1] + gamma_eff[1:]) * np.diff(times_array))
    reference = np.exp(-exposure)
    residual = np.asarray(edge_population) - reference
    return {
        "n_cells": n_cells,
        "n_modes": n,
        "profile": profile,
        "rates": rates.tolist(),
        "n_time": n_time,
        "times": times.tolist(),
        "gamma_eff": gamma_eff.tolist(),
        "integrated_exposure": exposure.tolist(),
        "edge_population": edge_population,
        "passive_reference": reference.tolist(),
        "residual": residual.tolist(),
        "no_jump_norm": no_jump_norm,
        "bulk_population": bulk_population,
        "max_abs_residual": float(np.max(np.abs(residual))),
        "max_bulk_population": float(np.max(np.abs(bulk_population))),
        "max_coherent_unitarity_residual": float(max(unitarity_norm)),
        "uniform_loss_analytic_control": profile == "uniform_loss",
        "interpretation": "Uniform loss is an exact scalar-loss control; edge-selective loss tests noncommuting loss/micromotion and open-system edge-bulk transfer.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", nargs="+", type=int, default=[2, 3, 4, 6])
    parser.add_argument("--profiles", nargs="+", choices=("edge_low_loss", "uniform_loss"), default=["edge_low_loss", "uniform_loss"])
    parser.add_argument("--n-time", type=int, default=257)
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "passive_single_particle_benchmark.json")
    args = parser.parse_args()
    if any(x < 1 for x in args.sizes) or args.n_time < 3:
        parser.error("sizes must be positive and n-time must be at least 3")
    results = [evaluate(n, profile, args.n_time) for n in args.sizes for profile in args.profiles]
    uniform = [r for r in results if r["profile"] == "uniform_loss"]
    output = {
        "status": "COMPLETED_DETERMINISTIC_PASSIVE_SINGLE_PARTICLE_BENCHMARK",
        "scientific_scope": "one-excitation passive dynamics; no Kerr, Fock truncation, trajectories or Wigner transform",
        "model": "same two-step chiral Floquet SSH backbone and coherent covariant edge vector",
        "reference": "n_edge(0) exp[- integral gamma_eff(t) dt]",
        "results": results,
        "uniform_loss_max_residual": float(max(r["max_abs_residual"] for r in uniform)),
        "edge_selective_max_residual": float(max(r["max_abs_residual"] for r in results if r["profile"] == "edge_low_loss")),
        "claim_boundary": [
            "isolates deterministic loss/micromotion effects",
            "does not validate Kerr dynamics or finite-Fock trajectories",
            "does not establish topological protection",
            "does not compute collective-mode Wigner negativity",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"uniform_loss_max_residual": output["uniform_loss_max_residual"], "edge_selective_max_residual": output["edge_selective_max_residual"], "rows": [{"n_cells": r["n_cells"], "profile": r["profile"], "max_abs_residual": r["max_abs_residual"], "max_bulk_population": r["max_bulk_population"]} for r in results]}, indent=2))


if __name__ == "__main__":
    main()
