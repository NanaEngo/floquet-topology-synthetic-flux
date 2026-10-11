#!/usr/bin/env python3
"""Audit paired edge-mode pilot outputs.

The independent unit is a matched disorder realization.  The script reports
sector contrasts for edge-mode population, local Wigner negativity and loss
currents.  It is descriptive and fail-closed: it never labels a protection gate
as passed from a single positive contrast.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

SECTORS = ("chiral_topological", "chiral_trivial")
PROFILES = ("edge_low_loss", "uniform_loss")
DISORDER_LEVELS = (0.0, 0.5, 1.0, 1.5)


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def mean_ci(values):
    x = np.asarray(values, dtype=float)
    n = len(x)
    mean = float(np.mean(x)) if n else float("nan")
    sd = float(np.std(x, ddof=1)) if n > 1 else 0.0
    half = 1.96 * sd / math.sqrt(n) if n else float("nan")
    return {
        "n": n, "mean": mean, "sd": sd,
        "ci95_low": mean - half, "ci95_high": mean + half,
        "positive_fraction": float(np.mean(x > 0)) if n else float("nan"),
    }


def audit(path: Path):
    data = json.loads(path.read_text())
    rows = data.get("records", [])
    n_mc = int(data["design"]["n_mc"])
    expected = 2 * len(PROFILES) * len(DISORDER_LEVELS) * n_mc
    if len(rows) != expected:
        raise ValueError(f"{path.name}: expected {expected} records, found {len(rows)}")
    if any(not finite(v) for row in rows for v in row.values() if isinstance(v, float)):
        raise ValueError(f"{path.name}: non-finite scalar detected")
    groups = {}
    for row in rows:
        key = (row["loss_profile"], float(row["disorder_level"]), int(row["sample"]), int(row["draw_seed"]))
        groups.setdefault(key, {})[row["sector"]] = row
    if any(set(g) != set(SECTORS) for g in groups.values()):
        raise ValueError(f"{path.name}: incomplete sector pairing")
    expected_groups = len(PROFILES) * len(DISORDER_LEVELS) * n_mc
    if len(groups) != expected_groups:
        raise ValueError(f"{path.name}: expected {expected_groups} paired units, found {len(groups)}")
    summaries = []
    for profile in PROFILES:
        for disorder in DISORDER_LEVELS:
            selected = [g for (p, d, _s, _seed), g in groups.items() if p == profile and d == disorder]
            metrics = {}
            for name in ("edge_mode_population", "negativity_site_0", "edge_loss_current", "bulk_loss_current"):
                values = [g["chiral_topological"][name] - g["chiral_trivial"][name] for g in selected]
                metrics[name + "_topological_minus_trivial"] = mean_ci(values)
            summaries.append({"loss_profile": profile, "disorder_level": disorder, **metrics})
    return {
        "input": str(path),
        "status": "STRUCTURAL_AUDIT_PASS",
        "scientific_status": "PROVISIONAL_EDGE_MODE_PILOT_AUDITED",
        "n_records": len(rows), "n_paired_units": len(groups),
        "design": data["design"], "summaries": summaries,
        "claim_boundary": [
            "descriptive finite-time matched dynamics only",
            "no topological-protection theorem",
            "no causal mechanism claim until passive bound and controls pass",
            "collective-mode Wigner negativity remains separately gated",
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    outputs = [audit(p) for p in args.input]
    result = {
        "status": "EDGE_MODE_PILOT_OUTPUTS_AUDITED",
        "files": outputs,
        "global_claim_boundary": "The pilot can establish implementation and effect-size diagnostics, not protection or mechanism by itself.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "files": len(outputs)}, indent=2))


if __name__ == "__main__":
    main()
