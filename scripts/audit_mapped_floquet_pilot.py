#!/usr/bin/env python3
"""Audit the mapped Floquet pilot and apply a fail-closed gate."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_mapped_floquet_dynamics import single_particle_matrix  # noqa: E402

PERIOD = 2.0 * np.pi
CANDIDATES = {"candidate_W1": (0.2, 1.2), "candidate_W0": (0.1, 0.1)}


def load(path):
    d = json.loads(path.read_text())
    return d, {(r["loss_profile"], float(r["disorder_level"]), r["sample"], r["sector"]): r
                for r in d["records"]}


def edge_diagnostic(n_cells, j1, j2):
    h0 = single_particle_matrix(n_cells, j1, j2, 0)
    h1 = single_particle_matrix(n_cells, j1, j2, 1)
    u = expm(-1j * h1 * np.pi) @ expm(-1j * h0 * np.pi)
    values, vectors = np.linalg.eig(u)
    edge_weights = []
    for i in range(len(values)):
        p = np.abs(vectors[:, i]) ** 2
        edge_weights.append(float(p[:2].sum() + p[-2:].sum()))
    return {
        "n_modes": 2 * n_cells,
        "max_edge_weight": float(max(edge_weights)),
        "mean_edge_weight": float(np.mean(edge_weights)),
        "n_modes_edge_weight_gt_0_5": int(sum(x > 0.5 for x in edge_weights)),
        "quasienergy_phases": sorted(float(x) for x in np.angle(values)),
        "interpretation": "Finite-chain edge weights are diagnostics only; they do not establish bulk-boundary correspondence."
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=Path, default=ROOT / "data" / "mapped_pilot")
    p.add_argument("--output", type=Path, default=ROOT / "data" / "mapped_pilot_analysis.json")
    args = p.parse_args()
    files = sorted(args.input_dir.glob("mapped_floquet_*cells.json"))
    if not files:
        raise SystemExit("no mapped pilot outputs found")
    summaries = []
    all_rows = {}
    for path in files:
        data, rows = load(path)
        cells = int(data["design"]["n_cells"])
        all_rows[cells] = rows
        for profile in ("edge_low_loss", "uniform_loss"):
            for disorder in (0.0, 0.5, 1.0, 1.5):
                diffs = []
                for sample in range(data["design"]["n_mc"]):
                    a = rows[(profile, disorder, sample, "candidate_W1")]
                    b = rows[(profile, disorder, sample, "candidate_W0")]
                    diffs.append(a["negativity_edge"] - b["negativity_edge"])
                summaries.append({
                    "n_cells": cells,
                    "loss_profile": profile,
                    "disorder_level": disorder,
                    "n_independent_realizations": len(diffs),
                    "W1_minus_W0_edge_negativity_mean": float(np.mean(diffs)),
                    "W1_minus_W0_edge_negativity_std": float(np.std(diffs, ddof=1)),
                })

    finite_size = []
    if 2 in all_rows and 3 in all_rows:
        keys = sorted(set(all_rows[2]) & set(all_rows[3]))
        for key in keys:
            finite_size.append(abs(all_rows[2][key]["negativity_edge"] - all_rows[3][key]["negativity_edge"]))
    edge = {
        str(cells): {
            name: edge_diagnostic(cells, *params)
            for name, params in CANDIDATES.items()
        }
        for cells in (2, 3, 4)
    }
    edge_profile = [x["W1_minus_W0_edge_negativity_mean"] for x in summaries
                    if x["loss_profile"] == "edge_low_loss"]
    output = {
        "status": "MAPPED_PILOT_AUDITED_FAIL_CLOSED",
        "inputs": [str(x) for x in files],
        "summaries": summaries,
        "single_particle_edge_diagnostics": edge,
        "finite_size": {
            "paired_records": len(finite_size),
            "max_abs_4mode_6mode_difference": max(finite_size) if finite_size else None,
            "mean_abs_4mode_6mode_difference": float(np.mean(finite_size)) if finite_size else None,
        },
        "gate": {
            "mapped_W1_outperforms_W0_in_all_edge_loss_conditions": bool(all(x > 0 for x in edge_profile)),
            "production_authorized": False,
            "reason": "The mapped pilot does not show the prespecified W1-positive finite-time contrast; the open-chain mapping and observable mechanism require revision before scaling n_mc."
        },
        "interpretation": "Pilot-only source-model evidence; no manuscript integration and no topological claim authorized.",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output["gate"], indent=2))
    print(json.dumps(output["finite_size"], indent=2))


if __name__ == "__main__":
    main()
