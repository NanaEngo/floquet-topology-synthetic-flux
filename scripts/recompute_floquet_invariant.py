#!/usr/bin/env python3
"""Recompute the framework Floquet winding and audit archived label mappings.

This is a provenance audit, not a proof of topology.  In particular, it checks
whether the two archived nonzero labels correspond to distinct parameter sets in
the available source code.  The result is intentionally kept separate from the
archived trace arrays.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from fl_qom.analysis.topological import TopologicalAnalyzer
from fl_qom.core.parameters import SystemParameters
from fl_qom.core.system import OptomechanicalSystem


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "trackA_invariant_audit.json"


def compute(j1: float, j2: float, mass: float = 1.0, delta: float = 0.5) -> int:
    system = OptomechanicalSystem(
        SystemParameters(omega_m=1.0e9), enable_logging=False, enable_metrics=False
    )
    analyzer = TopologicalAnalyzer(system, floquet_J1=j1, floquet_J2=j2)
    return int(analyzer.compute_floquet_winding_number(
        n_k=40, n_t=40, J1=j1, J2=j2, M=mass, delta=delta
    ))


def main() -> None:
    cases = {
        "archived_W_plus_half_label": {"J1": 0.15, "J2": 1.0, "declared_label": "+1/2"},
        "archived_W_minus_half_label": {"J1": 0.15, "J2": 1.0, "declared_label": "-1/2"},
        "archived_W_zero_label": {"J1": 0.8, "J2": 0.5, "declared_label": "0"},
        "n6_topological_source_control": {"J1": 0.1, "J2": 1.0, "declared_label": "topological"},
        "n6_trivial_source_control": {"J1": 1.0, "J2": 1.0, "declared_label": "trivial"},
    }
    results = {}
    for name, case in cases.items():
        results[name] = {
            **case,
            "recomputed_winding": compute(case["J1"], case["J2"]),
        }

    plus = results["archived_W_plus_half_label"]["recomputed_winding"]
    minus = results["archived_W_minus_half_label"]["recomputed_winding"]
    results_meta = {
        "status": "COMPLETED",
        "method": "FL_QOM TopologicalAnalyzer.compute_floquet_winding_number",
        "grid": {"n_k": 40, "n_t": 40},
        "software": {"numpy": np.__version__},
        "interpretation": (
            "The archived +/-1/2 labels use identical available J1/J2 parameters, "
            "so this audit cannot verify them as distinct invariant sectors."
        ),
        "archived_nonzero_labels_distinct_in_available_mapping": bool(plus != minus),
        "results": results,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(results_meta, indent=2) + "\n")
    print(json.dumps(results_meta, indent=2))


if __name__ == "__main__":
    main()
