#!/usr/bin/env python3
"""Audit Fock-cutoff convergence for the edge-mode pilot."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def finite(x):
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def audit(path: Path):
    d = json.loads(path.read_text())
    rows = d.get("records", [])
    expected = 2 * int(d["design"]["n_mc"])
    if len(rows) != expected:
        raise ValueError(f"{path.name}: expected {expected}, found {len(rows)}")
    groups = {}
    for r in rows:
        key = (r["sample"], r["draw_seed"])
        groups.setdefault(key, {})[r["sector"]] = r
        for name in ("edge_mode_population", "gamma_eff", "edge_loss_current", "bulk_loss_current"):
            if not finite(r[name]):
                raise ValueError(f"{path.name}: non-finite {name}")
    if any(set(g) != {"chiral_topological", "chiral_trivial"} for g in groups.values()):
        raise ValueError(f"{path.name}: broken sector pairing")
    deltas = [g["chiral_topological"]["edge_mode_population"] - g["chiral_trivial"]["edge_mode_population"] for g in groups.values()]
    passive = [r["max_abs_passive_residual"] for r in rows if r["passive_reference_status"] == "topological_passive_reference"]
    return {
        "file": str(path),
        "status": "PASS",
        "n_max": d["design"]["n_max"],
        "kerr": d["design"]["kerr"],
        "n_records": len(rows),
        "n_paired_units": len(groups),
        "edge_population_delta_mean": sum(deltas) / len(deltas),
        "edge_population_delta_min": min(deltas),
        "edge_population_delta_max": max(deltas),
        "max_passive_residual": (max(passive) if passive else None),
        "passive_reference_records": len(passive),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, nargs="+", required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    reports = [audit(x) for x in args.input]
    result = {
        "status": "EDGE_MODE_CONVERGENCE_AUDITED",
        "reports": reports,
        "interpretation": "Descriptive cutoff audit only; convergence does not establish physical calibration or topological protection.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "files": len(reports)}, indent=2))


if __name__ == "__main__":
    main()
