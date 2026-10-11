#!/usr/bin/env python3
"""Validate a standard chiral two-step Floquet SSH model.

The model is independent of the legacy FL_QOM winding routine.  Each drive step
contains only bipartite hopping, so Gamma=sigma_z chiral symmetry is exact.  The
script reports the two time-frame windings and their combinations
nu_0=(nu_A+nu_B)/2 and nu_pi=(nu_A-nu_B)/2, together with Floquet gaps and
finite-chain edge diagnostics.

This is a model-validation artifact, not yet a calculation of the archived
Wigner-negativity traces.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "chiral_floquet_phase_audit.json"
SIGMA_Z = np.diag([1.0, -1.0]).astype(complex)
PERIOD = 2.0 * np.pi
HALF_PERIOD = np.pi


def bulk_hamiltonian(k: float, intra: float, inter: float) -> np.ndarray:
    z = intra + inter * np.exp(-1j * k)
    return np.array([[0.0, z], [np.conjugate(z), 0.0]], dtype=complex)


def step_unitary(k: float, intra: float, inter: float, duration: float) -> np.ndarray:
    return expm(-1j * bulk_hamiltonian(k, intra, inter) * duration)


def physical_bulk_unitary(k: float, ja: float, jb: float) -> np.ndarray:
    ua = step_unitary(k, ja, 0.0, HALF_PERIOD)
    ub = step_unitary(k, 0.0, jb, HALF_PERIOD)
    return ub @ ua


def symmetric_frame_unitary(k: float, ja: float, jb: float, frame: str) -> np.ndarray:
    ah = step_unitary(k, ja, 0.0, HALF_PERIOD / 2.0)
    bh = step_unitary(k, 0.0, jb, HALF_PERIOD / 2.0)
    af = ah @ ah
    bf = bh @ bh
    if frame == "A":
        return bh @ af @ bh
    if frame == "B":
        return ah @ bf @ ah
    raise ValueError(frame)


def offdiag_winding(ja: float, jb: float, frame: str, n_k: int):
    ks = np.linspace(0.0, 2.0 * np.pi, n_k, endpoint=False)
    values = np.asarray([
        symmetric_frame_unitary(k, ja, jb, frame)[0, 1] for k in ks
    ])
    phases = np.unwrap(np.angle(values))
    winding_float = (phases[-1] - phases[0]) / (2.0 * np.pi)
    return {
        "winding_float": float(winding_float),
        "winding_integer": int(np.rint(winding_float)),
        "minimum_offdiagonal_magnitude": float(np.min(np.abs(values))),
    }


def floquet_gaps(ja: float, jb: float, n_k: int):
    gaps_zero = []
    gaps_pi = []
    chiral_residual_A = []
    chiral_residual_B = []
    for k in np.linspace(0.0, 2.0 * np.pi, n_k, endpoint=False):
        u = physical_bulk_unitary(k, ja, jb)
        phases = np.angle(np.linalg.eigvals(u))
        gaps_zero.extend(np.abs(phases))
        gaps_pi.extend(np.abs(np.abs(phases) - np.pi))
        ua = symmetric_frame_unitary(k, ja, jb, "A")
        ub = symmetric_frame_unitary(k, ja, jb, "B")
        chiral_residual_A.append(np.linalg.norm(SIGMA_Z @ ua @ SIGMA_Z - ua.conj().T))
        chiral_residual_B.append(np.linalg.norm(SIGMA_Z @ ub @ SIGMA_Z - ub.conj().T))
    return {
        "minimum_zero_gap": float(min(gaps_zero)),
        "minimum_pi_gap": float(min(gaps_pi)),
        "maximum_chiral_frame_A_residual": float(max(chiral_residual_A)),
        "maximum_chiral_frame_B_residual": float(max(chiral_residual_B)),
    }


def finite_chain_matrix(n_cells: int, intra: float, inter: float) -> np.ndarray:
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


def edge_diagnostic(n_cells: int, ja: float, jb: float, threshold: float = 0.5):
    ha = finite_chain_matrix(n_cells, ja, 0.0)
    hb = finite_chain_matrix(n_cells, 0.0, jb)
    u = expm(-1j * hb * HALF_PERIOD) @ expm(-1j * ha * HALF_PERIOD)
    values, vectors = np.linalg.eig(u)
    records = []
    for i, value in enumerate(values):
        p = np.abs(vectors[:, i]) ** 2
        edge_weight = float(p[:2].sum() + p[-2:].sum())
        phase = float(np.angle(value))
        distance_zero = abs(phase)
        distance_pi = abs(abs(phase) - np.pi)
        records.append({
            "phase": phase,
            "edge_weight": edge_weight,
            "near_zero": distance_zero < 0.15,
            "near_pi": distance_pi < 0.15,
        })
    return {
        "n_cells": n_cells,
        "threshold": threshold,
        "n_edge_localized": int(sum(x["edge_weight"] > threshold for x in records)),
        "n_edge_near_zero": int(sum(x["edge_weight"] > threshold and x["near_zero"] for x in records)),
        "n_edge_near_pi": int(sum(x["edge_weight"] > threshold and x["near_pi"] for x in records)),
        "max_edge_weight": float(max(x["edge_weight"] for x in records)),
        "states": records,
    }


def evaluate_point(ja: float, jb: float):
    rows = []
    for n in (40, 80, 120, 200):
        a = offdiag_winding(ja, jb, "A", n)
        b = offdiag_winding(ja, jb, "B", n)
        rows.append({
            "n_k": n,
            "frame_A": a,
            "frame_B": b,
            "nu_0": 0.5 * (a["winding_integer"] + b["winding_integer"]),
            "nu_pi": 0.5 * (a["winding_integer"] - b["winding_integer"]),
        })
    gaps = floquet_gaps(ja, jb, 200)
    stable = all(x["nu_0"] == rows[0]["nu_0"] and x["nu_pi"] == rows[0]["nu_pi"] for x in rows)
    return {
        "ja_intra_step_A": ja,
        "jb_inter_step_B": jb,
        "resolution_results": rows,
        "floquet_gaps": gaps,
        "resolution_stable": bool(stable),
        "chiral_unitary_validated": bool(gaps["maximum_chiral_frame_A_residual"] < 1e-10 and gaps["maximum_chiral_frame_B_residual"] < 1e-10),
        "bulk_point_candidate_pass": bool(stable and gaps["minimum_zero_gap"] > 1e-8 and gaps["minimum_pi_gap"] > 1e-8 and gaps["maximum_chiral_frame_A_residual"] < 1e-10 and gaps["maximum_chiral_frame_B_residual"] < 1e-10),
        "edge_diagnostics": {str(n): edge_diagnostic(n, ja, jb) for n in (4, 8, 12)},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--n-points", type=int, default=15)
    args = parser.parse_args()
    values = np.linspace(0.1, 1.5, args.n_points)
    points = []
    for ja in values:
        for jb in values:
            points.append(evaluate_point(float(ja), float(jb)))
    selected = {
        "candidate_topological": evaluate_point(0.1, 1.2),
        "candidate_trivial": evaluate_point(1.2, 0.1),
    }
    output = {
        "status": "COMPLETED_INDEPENDENT_CHIRAL_FLOQUET_AUDIT",
        "model": {
            "type": "two-step chiral Floquet SSH",
            "period": PERIOD,
            "step_A": "pure intracell hopping ja for duration pi",
            "step_B": "pure intercell hopping jb for duration pi",
            "chiral_operator": "sigma_z",
        },
        "invariant_definition": {
            "frame_windings": "w_A and w_B from the phase winding of the off-diagonal block of symmetric time-frame unitaries",
            "nu_0": "(w_A+w_B)/2",
            "nu_pi": "(w_A-w_B)/2",
        },
        "selected_points": selected,
        "scan": {
            "n_points_per_axis": args.n_points,
            "J_range": [float(values[0]), float(values[-1])],
            "n_records": len(points),
            "stable_records": sum(x["bulk_point_candidate_pass"] for x in points),
            "records": points,
        },
        "interpretation": "This independent chiral model provides a mathematically controlled topology test. It is not yet connected to the archived Wigner-negativity arrays or to the previous static-chain simulations.",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({
        "selected_topological": {k: selected["candidate_topological"][k] for k in ["resolution_stable", "floquet_gaps", "bulk_point_candidate_pass"]},
        "selected_trivial": {k: selected["candidate_trivial"][k] for k in ["resolution_stable", "floquet_gaps", "bulk_point_candidate_pass"]},
        "stable_scan_records": output["scan"]["stable_records"],
    }, indent=2))


if __name__ == "__main__":
    main()
