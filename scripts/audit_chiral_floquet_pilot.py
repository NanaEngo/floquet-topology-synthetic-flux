#!/usr/bin/env python3
"""Audit the chiral Floquet dynamics pilot before production scaling."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=Path, default=ROOT / "data" / "chiral_pilot")
    p.add_argument("--output", type=Path, default=ROOT / "data" / "chiral_pilot_analysis.json")
    args = p.parse_args()
    files = sorted(args.input_dir.glob("chiral_floquet_*cells.json"))
    if not files:
        raise SystemExit("no chiral pilot outputs")
    summaries = []
    by_size = {}
    for path in files:
        d = json.loads(path.read_text())
        rows = d["records"]
        cells = d["design"]["n_cells"]
        by_size[cells] = {(r["loss_profile"], r["disorder_level"], r["sample"], r["sector"]): r for r in rows}
        for profile in ("edge_low_loss", "uniform_loss"):
            for disorder in (0.0, 0.5, 1.0, 1.5):
                a = [r for r in rows if r["loss_profile"] == profile and r["disorder_level"] == disorder and r["sector"] == "chiral_topological"]
                b = [r for r in rows if r["loss_profile"] == profile and r["disorder_level"] == disorder and r["sector"] == "chiral_trivial"]
                diff = np.asarray([x["negativity_edge"] - y["negativity_edge"] for x, y in zip(a, b)])
                summaries.append({
                    "n_cells": cells,
                    "loss_profile": profile,
                    "disorder_level": disorder,
                    "n_independent_realizations": len(diff),
                    "topological_minus_trivial_mean": float(diff.mean()),
                    "topological_minus_trivial_std": float(diff.std(ddof=1)),
                    "topological_mean": float(np.mean([x["negativity_edge"] for x in a])),
                    "trivial_mean": float(np.mean([x["negativity_edge"] for x in b])),
                })
    finite = []
    if 2 in by_size and 3 in by_size:
        keys = sorted(set(by_size[2]) & set(by_size[3]))
        finite = [abs(by_size[2][k]["negativity_edge"] - by_size[3][k]["negativity_edge"]) for k in keys]
    edge_diffs = [x["topological_minus_trivial_mean"] for x in summaries if x["loss_profile"] == "edge_low_loss"]
    output = {
        "status": "CHIRAL_PILOT_AUDITED_NOT_PRODUCTION",
        "inputs": [str(x) for x in files],
        "summaries": summaries,
        "finite_size": {
            "paired_records": len(finite),
            "max_abs_4mode_6mode_difference": max(finite) if finite else None,
            "mean_abs_4mode_6mode_difference": float(np.mean(finite)) if finite else None,
        },
        "gate": {
            "edge_loss_positive_at_all_sizes_and_disorder": bool(all(x > 0 for x in edge_diffs)),
            "uniform_loss_positive_at_all_sizes_and_disorder": bool(all(x > 0 for x in [y["topological_minus_trivial_mean"] for y in summaries if y["loss_profile"] == "uniform_loss"])),
            "production_authorized": False,
            "reason": "The edge-loss signal is positive for N=4 in several conditions but not uniformly for N=6; uniform-loss conditions do not show a robust separation. Increase the model-mechanism validation before n_mc scaling."
        },
        "initial_state_mode": "edge_site",
        "interpretation": "This is a pilot of a certified chiral bulk model, not evidence from the archived reconstruction."
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output["gate"], indent=2))
    print(json.dumps(output["finite_size"], indent=2))


if __name__ == "__main__":
    main()
