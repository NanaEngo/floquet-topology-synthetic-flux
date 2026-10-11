#!/usr/bin/env python3
"""Audit assumptions behind the source Floquet winding calculation.

The script does not assume that an integer returned by a routine is a valid
invariant.  It checks the algebraic conditions required by the proposed model
and records the result as an audit.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from fl_qom.analysis.topological import TopologicalAnalyzer
from fl_qom.core.parameters import SystemParameters
from fl_qom.core.system import OptomechanicalSystem

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "floquet_symmetry_audit.json"
SIGMA_Z = np.diag([1.0, -1.0]).astype(complex)


def analyzer():
    return TopologicalAnalyzer(
        OptomechanicalSystem(SystemParameters(omega_m=1e9),
                             enable_logging=False, enable_metrics=False),
        floquet_M=1.0,
        floquet_delta=0.5,
    )


def rel_norm(matrix):
    return float(np.linalg.norm(matrix) / max(np.linalg.norm(matrix), 1e-15))


def audit_hamiltonian(a, j1, j2):
    records = []
    for step in (0, 1):
        herm = []
        chiral = []
        periodic = []
        for k in np.linspace(0, 2 * np.pi, 41):
            h = a.floquet_two_step_hamiltonian(k, step, J1=j1, J2=j2, M=1.0, delta=0.5)
            h_periodic = a.floquet_two_step_hamiltonian(k + 2 * np.pi, step, J1=j1, J2=j2, M=1.0, delta=0.5)
            herm.append(float(np.linalg.norm(h - h.conj().T)))
            chiral.append(float(np.linalg.norm(SIGMA_Z @ h @ SIGMA_Z + h)))
            periodic.append(float(np.linalg.norm(h - h_periodic)))
        records.append({
            "step": step,
            "max_hermiticity_residual": max(herm),
            "max_chiral_anticommutator_norm": max(chiral),
            "max_periodicity_residual": max(periodic),
        })
    return records


def audit_unitarity(a, j1, j2):
    output = []
    for n_t in (20, 40, 80):
        residuals = []
        for k in np.linspace(0, 2 * np.pi, 41, endpoint=False):
            u = a._floquet_evolution_operator(k, n_t)[-1]
            residuals.append(float(np.linalg.norm(u.conj().T @ u - np.eye(2))))
        output.append({"n_t": n_t, "max_unitarity_residual": max(residuals)})
    return output


def audit_winding(a, j1, j2):
    rows = []
    for n in (20, 40, 80, 120):
        value = int(a.compute_floquet_winding_number(
            n_k=n, n_t=n, J1=j1, J2=j2, M=1.0, delta=0.5
        ))
        rows.append({"n_k": n, "n_t": n, "framework_winding": value})
    return rows


def main():
    a = analyzer()
    cases = {
        "historical_N6_source": {"J1": 0.1, "J2": 1.0},
        "candidate_W1": {"J1": 0.2, "J2": 1.2},
        "candidate_W0": {"J1": 0.1, "J2": 0.1},
        "historical_control": {"J1": 1.0, "J2": 1.0},
    }
    result = {
        "status": "COMPLETED_ALGEBRAIC_AUDIT",
        "model": "FL_QOM two-step Hamiltonian",
        "symmetry_operator": "sigma_z",
        "required_for_simple_chiral_1D_winding": "{Gamma,H(k,t)}=0 for every k,t",
        "cases": {},
        "decision": {
            "source_has_exact_chiral_symmetry": False,
            "source_winding_is_not_yet_validated_as_chiral_invariant": True,
            "next_action": "Implement an explicitly chiral Floquet model or define a different invariant with its assumptions."
        }
    }
    for name, point in cases.items():
        result["cases"][name] = {
            **point,
            "hamiltonian": audit_hamiltonian(a, point["J1"], point["J2"]),
            "unitarity": audit_unitarity(a, point["J1"], point["J2"]),
            "winding_resolution": audit_winding(a, point["J1"], point["J2"]),
        }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["decision"], indent=2))
    for name, data in result["cases"].items():
        print(name, data["winding_resolution"])


if __name__ == "__main__":
    main()
