#!/usr/bin/env python3
"""Scan the available two-step FL_QOM model before selecting final sectors.

This script is an audit of the *implemented* source model.  It is deliberately
not named or presented as a proof of a topological phase diagram.  For each
(J1, J2) point it records the framework winding output at two resolutions and a
separate period-operator quasienergy-gap diagnostic.  A point is marked stable
only when the two reported winding resolutions agree and the gap is non-zero.

The scan is written to the V2608 workspace and never modifies archived arrays.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from fl_qom.analysis.topological import TopologicalAnalyzer
from fl_qom.core.parameters import SystemParameters
from fl_qom.core.system import OptomechanicalSystem


DEFAULT_J1 = np.linspace(0.1, 1.5, 15)
DEFAULT_J2 = np.linspace(0.1, 1.5, 15)
MASS = 1.0
DELTA = 0.5


def make_analyzer() -> TopologicalAnalyzer:
    system = OptomechanicalSystem(
        SystemParameters(omega_m=1.0e9),
        enable_logging=False,
        enable_metrics=False,
    )
    return TopologicalAnalyzer(system, floquet_M=MASS, floquet_delta=DELTA)


def period_operator_at_k(k: float, j1: float, j2: float, n_t: int) -> np.ndarray:
    analyzer = make_analyzer()
    dt = 2.0 * np.pi / n_t
    h0 = analyzer.floquet_two_step_hamiltonian(
        k, 0, J1=j1, J2=j2, M=MASS, delta=DELTA
    )
    h1 = analyzer.floquet_two_step_hamiltonian(
        k, 1, J1=j1, J2=j2, M=MASS, delta=DELTA
    )
    u = np.eye(2, dtype=complex)
    for m in range(n_t):
        h = h0 if m % 2 == 0 else h1
        u = expm(-1j * h * dt) @ u
    return u


def quasienergy_gap(j1: float, j2: float, n_k: int, n_t: int) -> float:
    """Minimum circular eigenphase separation of the period operator."""
    gaps = []
    for k in np.linspace(0.0, 2.0 * np.pi, n_k, endpoint=False):
        evals = np.linalg.eigvals(period_operator_at_k(k, j1, j2, n_t))
        phases = np.sort(np.angle(evals))
        separation = abs(float(phases[1] - phases[0]))
        gaps.append(min(separation, 2.0 * np.pi - separation))
    return float(min(gaps))


def evaluate_point(j1: float, j2: float) -> dict:
    analyzer = make_analyzer()
    w20 = int(analyzer.compute_floquet_winding_number(
        n_k=20, n_t=20, J1=j1, J2=j2, M=MASS, delta=DELTA
    ))
    w40 = int(analyzer.compute_floquet_winding_number(
        n_k=40, n_t=40, J1=j1, J2=j2, M=MASS, delta=DELTA
    ))
    gap20 = quasienergy_gap(j1, j2, n_k=20, n_t=20)
    gap40 = quasienergy_gap(j1, j2, n_k=40, n_t=40)
    return {
        "J1": float(j1),
        "J2": float(j2),
        "M": MASS,
        "delta": DELTA,
        "framework_winding_n20": w20,
        "framework_winding_n40": w40,
        "resolution_agreement": bool(w20 == w40),
        "min_quasienergy_gap_n20": gap20,
        "min_quasienergy_gap_n40": gap40,
        "gap_positive_at_both_resolutions": bool(gap20 > 1e-8 and gap40 > 1e-8),
        "audit_stable_label": bool(w20 == w40 and gap20 > 1e-8 and gap40 > 1e-8),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1]
                        / "data" / "floquet_phase_audit.json")
    parser.add_argument("--n-points", type=int, default=15)
    args = parser.parse_args()
    if args.n_points < 3:
        parser.error("--n-points must be at least 3")

    j1_values = np.linspace(DEFAULT_J1[0], DEFAULT_J1[-1], args.n_points)
    j2_values = np.linspace(DEFAULT_J2[0], DEFAULT_J2[-1], args.n_points)
    records = [evaluate_point(float(j1), float(j2))
               for j1 in j1_values for j2 in j2_values]

    stable = [r for r in records if r["audit_stable_label"]]
    winding_values = sorted({r["framework_winding_n40"] for r in records})
    output = {
        "status": "COMPLETED_AUDIT_SCAN",
        "interpretation": (
            "This is a source-model audit. Agreement of framework winding values "
            "and a positive quasienergy gap are necessary diagnostics, not an "
            "independent proof of an open-system topological invariant."
        ),
        "model": {
            "implementation": "FL_QOM TopologicalAnalyzer.floquet_two_step_hamiltonian",
            "M": MASS,
            "delta": DELTA,
            "j1_range": [float(j1_values[0]), float(j1_values[-1])],
            "j2_range": [float(j2_values[0]), float(j2_values[-1])],
            "n_points_per_axis": args.n_points,
        },
        "diagnostics": {
            "n_records": len(records),
            "stable_label_records": len(stable),
            "framework_winding_values_n40": winding_values,
            "stable_records_by_winding": {
                str(w): sum(r["audit_stable_label"] and r["framework_winding_n40"] == w
                            for r in records)
                for w in winding_values
            },
        },
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output["diagnostics"], indent=2))
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
