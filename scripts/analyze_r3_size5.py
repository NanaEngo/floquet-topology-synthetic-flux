#!/usr/bin/env python3
"""Score the R3 size-5 run against the frozen out-of-sample prediction.

Implements the L3 gate of AUDIT_RESEARCH_REFINEMENTS_20261001.md (line 155):
the 5-cell (10-mode) edge negativity contrast measured by the R3 run is
compared against ``data/r3_size5_prediction.json``, which was frozen before
the run was submitted.  The prediction file is never rewritten by this
script.

Outputs:
    data/r3_size5_analysis.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

R3_FILE = ROOT / "data" / "chiral_production" / "chiral_floquet_5cells_NMC300.json"
PRED_FILE = ROOT / "data" / "r3_size5_prediction.json"
PROD_SUMMARY = ROOT / "data" / "chiral_production_summary.csv"
OBS = "negativity_edge_topological_minus_trivial"


def paired_contrast(records: list, disorder: float, profile: str):
    """Paired topological-minus-trivial edge negativity contrast per sample.

    Args:
        records (list): R3 record rows.
        disorder (float): Disorder level filter (dimensionless).
        profile (str): Loss profile filter.

    Returns:
        tuple: (samples, contrasts) arrays (dimensionless).
    """
    sel = [r for r in records if r["disorder_level"] == disorder
           and r["loss_profile"] == profile]
    topo = {r["sample"]: r["negativity_edge"] for r in sel
            if r["sector"] == "chiral_topological"}
    triv = {r["sample"]: r["negativity_edge"] for r in sel
            if r["sector"] == "chiral_trivial"}
    common = sorted(set(topo) & set(triv))
    return (np.array(common),
            np.array([topo[s] - triv[s] for s in common]))


def production_mean(disorder: float, profile: str) -> dict:
    """Production contrast statistics from the summary CSV.

    Args:
        disorder (float): Disorder level (dimensionless).
        profile (str): Loss profile.

    Returns:
        dict: ``{n_cells: {"mean": ..., "se": ...}}``.
    """
    import csv
    out = {}
    with PROD_SUMMARY.open() as fh:
        for row in csv.DictReader(fh):
            if (row["observable"] == OBS and float(row["disorder_level"]) == disorder
                    and row["loss_profile"] == profile):
                out[int(row["n_cells"])] = {"mean": float(row["mean"]),
                                            "se": float(row["se"])}
    return out


def main() -> None:
    """Aggregate the R3 contrasts and score them against the frozen forecast."""
    r3 = json.loads(R3_FILE.read_text())
    pred = json.loads(PRED_FILE.read_text())
    records = r3["records"]

    results = {}
    for disorder in sorted({r["disorder_level"] for r in records}):
        for profile in sorted({r["loss_profile"] for r in records}):
            samples, d = paired_contrast(records, disorder, profile)
            mean = float(d.mean())
            se = float(d.std(ddof=1) / np.sqrt(len(d)))
            key = f"disorder_{disorder}_{profile}"
            frozen = pred["predictions"][key]["ensemble_prediction"]
            measured = {
                "n": int(len(d)), "mean": mean, "se": se,
                "ci95": [mean - 1.96 * se, mean + 1.96 * se],
                "positive_fraction": float(np.mean(d > 0)),
            }
            sign = ("positive" if mean - 1.96 * se > 0 else
                    "negative" if mean + 1.96 * se < 0 else "unresolved")
            hit_interval = frozen["ci95_low"] <= mean <= frozen["ci95_high"]
            hit_sign = sign == pred["predictions"][key]["sign_prediction"]
            results[key] = {
                "measured": measured, "measured_sign": sign,
                "frozen_prediction": frozen,
                "frozen_sign": pred["predictions"][key]["sign_prediction"],
                "inside_frozen_interval": bool(hit_interval),
                "sign_match": bool(hit_sign),
                "model_predictions": {m: v["pred"] for m, v in
                                      pred["predictions"][key]["models"].items()},
            }

    # size trend table at d = 0, edge_low (headline condition)
    trend = production_mean(0.0, "edge_low_loss")
    s5, d5 = paired_contrast(records, 0.0, "edge_low_loss")
    trend[5] = {"mean": float(d5.mean()),
                "se": float(d5.std(ddof=1) / np.sqrt(len(d5)))}

    out = {
        "status": "COMPLETED_R3_SIZE5_ANALYSIS",
        "audit_source": "AUDIT_RESEARCH_REFINEMENTS_20261001.md (L3, line 155)",
        "prediction_file": str(PRED_FILE.relative_to(ROOT)),
        "prediction_frozen_unix": pred.get("created_unix"),
        "scoring": {
            "all_inside_interval": all(v["inside_frozen_interval"] for v in results.values()),
            "all_sign_match": all(v["sign_match"] for v in results.values()),
        },
        "conditions": results,
        "size_trend_d0_edge_low": {str(k): v for k, v in sorted(trend.items())},
    }
    dest = ROOT / "data" / "r3_size5_analysis.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dest}")
    for k, v in results.items():
        m = v["measured"]
        print(f"  {k}: {m['mean']:+.6f} ± {m['se']:.6f} "
              f"[{m['ci95'][0]:+.6f},{m['ci95'][1]:+.6f}] -> {v['measured_sign']} "
              f"(frozen: {v['frozen_prediction']['mean']:+.6f} "
              f"[{v['frozen_prediction']['ci95_low']:+.6f},"
              f"{v['frozen_prediction']['ci95_high']:+.6f}] -> {v['frozen_sign']}) "
              f"interval_hit={v['inside_frozen_interval']} sign_match={v['sign_match']}")
    print("size trend (d=0, edge_low):",
          {k: round(v["mean"], 6) for k, v in sorted(trend.items())})


if __name__ == "__main__":
    main()
