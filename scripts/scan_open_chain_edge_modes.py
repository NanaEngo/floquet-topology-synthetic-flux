#!/usr/bin/env python3
"""Audit open-chain edge localization for the source-model candidate scan.

This deliberately tests the proposed bulk-to-edge interpretation independently
of Wigner negativity.  It is a finite-chain diagnostic, not a proof of
bulk-boundary correspondence.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from run_mapped_floquet_dynamics import single_particle_matrix

ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "data" / "floquet_phase_audit.json"
OUTPUT = ROOT / "data" / "open_chain_edge_audit.json"
CELL_SIZES = (4, 8, 12)


def audit_point(cells, j1, j2):
    h0 = single_particle_matrix(cells, j1, j2, 0)
    h1 = single_particle_matrix(cells, j1, j2, 1)
    u = expm(-1j * h1 * np.pi) @ expm(-1j * h0 * np.pi)
    eigenvalues, vectors = np.linalg.eig(u)
    weights = []
    for i in range(len(eigenvalues)):
        p = np.abs(vectors[:, i]) ** 2
        weights.append(float(p[:2].sum() + p[-2:].sum()))
    return {
        "max_edge_weight": float(max(weights)),
        "mean_edge_weight": float(np.mean(weights)),
        "n_modes_edge_weight_gt_0_5": int(sum(x > 0.5 for x in weights)),
        "quasienergy_phases": sorted(float(x) for x in np.angle(eigenvalues)),
    }


def main():
    scan = json.loads(SCAN.read_text())
    rows = [r for r in scan["records"] if r["audit_stable_label"]]
    records = []
    for r in rows:
        for cells in CELL_SIZES:
            records.append({
                "J1": r["J1"], "J2": r["J2"],
                "framework_winding": r["framework_winding_n40"],
                "n_cells": cells,
                **audit_point(cells, r["J1"], r["J2"]),
            })
    summary = {}
    for cells in CELL_SIZES:
        summary[str(cells)] = {}
        for winding in (0, 1):
            sub = [r for r in records if r["n_cells"] == cells and r["framework_winding"] == winding]
            summary[str(cells)][str(winding)] = {
                "n_points": len(sub),
                "max_of_max_edge_weight": max(r["max_edge_weight"] for r in sub),
                "mean_max_edge_weight": float(np.mean([r["max_edge_weight"] for r in sub])),
                "points_with_gt_0_5_edge_modes": sum(r["n_modes_edge_weight_gt_0_5"] > 0 for r in sub),
            }
    output = {
        "status": "COMPLETED_EDGE_DIAGNOSTIC_NOT_TOPOLOGY_PROOF",
        "input_scan": "data/floquet_phase_audit.json",
        "cell_sizes": CELL_SIZES,
        "edge_definition": "probability on the first and last single-particle sites",
        "threshold": 0.5,
        "summary": summary,
        "records": records,
        "interpretation": "The diagnostic tests, but does not establish, bulk-boundary correspondence. The W=1 source-model label does not show a uniformly stronger edge signature than W=0 across finite sizes.",
        "gate": {
            "mechanism_supported": False,
            "production_authorized": False,
            "reason": "No size-robust edge-localization separation was established for the selected mapping."
        }
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"summary": summary, "gate": output["gate"]}, indent=2))


if __name__ == "__main__":
    main()
