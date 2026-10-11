#!/usr/bin/env python3
"""Validate provisional source-model sector candidates across resolutions.

The output is intentionally a candidate audit.  It does not certify a physical
open-system topological invariant or establish correspondence with the QuTiP
N-site dynamics model.
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
OUTPUT = ROOT / "data" / "candidate_sector_validation.json"
MASS = 1.0
DELTA = 0.5
CANDIDATES = {
    "candidate_W1": {"J1": 0.2, "J2": 1.2, "provisional_label": 1},
    "candidate_W0": {"J1": 0.1, "J2": 0.1, "provisional_label": 0},
}
RESOLUTIONS = (20, 40, 80, 120)


def analyzer() -> TopologicalAnalyzer:
    system = OptomechanicalSystem(
        SystemParameters(omega_m=1.0e9),
        enable_logging=False,
        enable_metrics=False,
    )
    return TopologicalAnalyzer(system, floquet_M=MASS, floquet_delta=DELTA)


def gap_at_resolution(j1: float, j2: float, n: int) -> float:
    a = analyzer()
    dt = 2.0 * np.pi / n
    gaps = []
    for k in np.linspace(0.0, 2.0 * np.pi, n, endpoint=False):
        h0 = a.floquet_two_step_hamiltonian(k, 0, J1=j1, J2=j2, M=MASS, delta=DELTA)
        h1 = a.floquet_two_step_hamiltonian(k, 1, J1=j1, J2=j2, M=MASS, delta=DELTA)
        u = np.eye(2, dtype=complex)
        for step in range(n):
            u = expm(-1j * (h0 if step % 2 == 0 else h1) * dt) @ u
        phases = np.sort(np.angle(np.linalg.eigvals(u)))
        sep = abs(float(phases[1] - phases[0]))
        gaps.append(min(sep, 2.0 * np.pi - sep))
    return float(min(gaps))


def main() -> None:
    output = {
        "status": "CANDIDATES_AUDITED_NOT_CERTIFIED",
        "model_manifest": "data/model_manifest.json",
        "resolutions": list(RESOLUTIONS),
        "candidates": {},
        "decision_rule": "Candidate is source-model-stable only if all winding outputs agree and all minimum gaps are >1e-8.",
    }
    for name, point in CANDIDATES.items():
        values = []
        for n in RESOLUTIONS:
            a = analyzer()
            winding = int(a.compute_floquet_winding_number(
                n_k=n, n_t=n, J1=point["J1"], J2=point["J2"], M=MASS, delta=DELTA
            ))
            values.append({
                "resolution": n,
                "framework_winding": winding,
                "min_quasienergy_gap": gap_at_resolution(point["J1"], point["J2"], n),
            })
        windings = {x["framework_winding"] for x in values}
        gaps = [x["min_quasienergy_gap"] for x in values]
        output["candidates"][name] = {
            **point,
            "resolution_results": values,
            "winding_values": sorted(windings),
            "resolution_stable": len(windings) == 1,
            "all_gaps_positive": min(gaps) > 1e-8,
            "source_model_candidate_pass": len(windings) == 1 and min(gaps) > 1e-8,
            "physical_sector_certified": False,
            "reason_not_certified": "Mapping to the open dissipative N-site dynamics remains unresolved.",
        }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
