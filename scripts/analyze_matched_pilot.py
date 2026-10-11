#!/usr/bin/env python3
"""Analyze the paired N=4/N=6 pilot and expose design failures early."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load(path: Path):
    data = json.loads(path.read_text())
    by_key = {}
    for row in data["records"]:
        key = (row["loss_profile"], float(row["disorder_level"]), row["sample"])
        by_key[(key, row["sector"])] = row
    return data, by_key


def bootstrap_ci(values, seed=20260814, n_boot=5000):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    means = np.mean(rng.choice(values, size=(n_boot, len(values)), replace=True), axis=1)
    return [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=ROOT / "data" / "pilot")
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "pilot_analysis.json")
    args = parser.parse_args()

    files = sorted(args.input_dir.glob("matched_dynamics_N*_pilot.json"))
    if not files:
        raise SystemExit("no pilot outputs found")
    summaries = []
    raw_by_size = {}
    for path in files:
        data, rows = load(path)
        n_sites = data["design"]["n_sites"]
        raw_by_size[n_sites] = rows
        for profile in ("edge_low_loss", "uniform_loss"):
            for disorder in (0.0, 0.5, 1.0, 1.5):
                diffs = []
                edge_w1 = []
                center_w1 = []
                for sample in range(data["design"]["n_mc"]):
                    a = rows[((profile, disorder, sample), "candidate_W1")]
                    b = rows[((profile, disorder, sample), "candidate_W0")]
                    diffs.append(a["negativity_edge"] - b["negativity_edge"])
                    edge_w1.append(a["edge_population"])
                    center_w1.append(a["negativity_center"])
                summaries.append({
                    "n_sites": n_sites,
                    "loss_profile": profile,
                    "disorder_level": disorder,
                    "n_independent_realizations": len(diffs),
                    "paired_edge_negativity_difference_mean": float(np.mean(diffs)),
                    "paired_edge_negativity_difference_std": float(np.std(diffs, ddof=1)),
                    "paired_edge_negativity_difference_bootstrap95": bootstrap_ci(diffs),
                    "W1_edge_population_mean": float(np.mean(edge_w1)),
                    "W1_center_negativity_mean": float(np.mean(center_w1)),
                })

    finite_size = []
    if 4 in raw_by_size and 6 in raw_by_size:
        keys = sorted(set(k for k in raw_by_size[4] if k in raw_by_size[6]))
        for key in keys:
            r4 = raw_by_size[4][key]
            r6 = raw_by_size[6][key]
            finite_size.append({
                "key": list(key),
                "sector": key[1] if len(key) > 1 else "unknown",
                "edge_negativity_abs_difference": abs(r4["negativity_edge"] - r6["negativity_edge"]),
                "edge_population_abs_difference": abs(r4["edge_population"] - r6["edge_population"]),
            })
    max_size_difference = max((x["edge_negativity_abs_difference"] for x in finite_size), default=None)
    output = {
        "status": "PILOT_ANALYZED_NOT_MANUSCRIPT_RESULT",
        "inputs": [str(p) for p in files],
        "summaries": summaries,
        "finite_size_comparison": {
            "paired_records": len(finite_size),
            "max_abs_N4_N6_edge_negativity_difference": max_size_difference,
            "interpretation": "Near-identical N=4/N=6 outputs require a longer-time or less strongly lossy finite-size design before claiming size dependence."
        },
        "alerts": [
            "The W1-W0 contrast changes sign across loss/disorder conditions in the pilot.",
            "The N=4 and N=6 outputs are nearly identical over the current t<=1 protocol.",
            "The pilot does not validate the mapping between FL_QOM candidate labels and the open-chain dynamics."
        ],
        "decision": "DO_NOT_LAUNCH_N_MC_500_UNCHANGED",
        "recommended_revision": "Add a prespecified longer-time or parameter-matched transport diagnostic, then repeat a small pilot before production."
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output["finite_size_comparison"], indent=2))
    print(output["decision"])


if __name__ == "__main__":
    main()
