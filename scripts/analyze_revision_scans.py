#!/usr/bin/env python3
"""Aggregate the 2026-08-17 revision scans into a single machine-readable summary.

Reads the per-scan outputs in data/revision_scans/ (merged JSON or, during
execution, the per-(profile,disorder) fragment files) and assembles matched
paired contrasts keyed by (scan, size, cutoff/ratio/period, profile, disorder).
Contrast = chiral_topological - chiral_trivial within a matched disorder draw,
then summarized over the independent units.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCAN_DIR = ROOT / "data" / "revision_scans"
SECTORS = ("chiral_topological", "chiral_trivial")


def load_records(path: Path) -> list[dict]:
    if path.is_file():
        data = json.loads(path.read_text())
        if isinstance(data, dict) and "records" in data:
            return data["records"]
        return data
    fdir = path.with_suffix("") / "_fragments"
    rows = []
    for frag in sorted(fdir.glob("*.json")):
        rows += json.loads(frag.read_text())
    if rows:
        return rows
    raise FileNotFoundError(f"no records for {path}")


def summarize(rows: list[dict], n_cells: int, profile: str, disorder: float):
    groups = {}
    for r in rows:
        key = (r["loss_profile"], float(r["disorder_level"]), int(r["sample"]))
        groups.setdefault(key, {})[r["sector"]] = r
    vals = []
    for (p, d, _), pair in groups.items():
        if p == profile and d == disorder:
            if set(pair) != set(SECTORS):
                raise ValueError("incomplete sector pairing")
            vals.append(pair["chiral_topological"]["negativity_edge"]
                        - pair["chiral_trivial"]["negativity_edge"])
    vals = np.asarray(vals)
    return {
        "n": int(len(vals)),
        "mean": float(vals.mean()),
        "sd": float(vals.std(ddof=1)) if len(vals) > 1 else 0.0,
        "positive_fraction": float((vals > 0).mean()) if len(vals) else float("nan"),
    }


def stats(values):
    """Mean, sd, standard error, 95% normal interval and sign fraction."""
    x = np.asarray(list(values), dtype=float)
    n = int(x.size)
    if n == 0:
        return {"n": 0, "mean": float("nan"), "sd": float("nan"), "se": float("nan"),
                "ci95_low": float("nan"), "ci95_high": float("nan"),
                "positive_fraction": float("nan")}
    mean = float(x.mean())
    sd = float(x.std(ddof=1)) if n > 1 else 0.0
    se = sd / np.sqrt(n) if n else float("nan")
    return {
        "n": n, "mean": mean, "sd": sd, "se": float(se),
        "ci95_low": float(mean - 1.96 * se), "ci95_high": float(mean + 1.96 * se),
        "positive_fraction": float((x > 0).mean()),
    }


def condition_groups(rows: list[dict]):
    """Group records into matched sector pairs keyed by (profile, disorder)."""
    groups = {}
    for r in rows:
        key = (r["loss_profile"], float(r["disorder_level"]), int(r["sample"]))
        groups.setdefault(key, {})[r["sector"]] = r
    return groups


def summarize_micromotion(path: Path) -> dict:
    """Contrast-versus-time curves for the sub-period (R7) scan."""
    rows = load_records(path)
    groups = condition_groups(rows)
    out = {}
    for profile, disorder in sorted({(p, d) for (p, d, _ ) in groups}):
        pairs = [g for (p, d, _), g in groups.items()
                 if p == profile and d == disorder and set(g) == set(SECTORS)]
        if not pairs:
            continue
        n_t = len(next(iter(pairs))["chiral_topological"]["time_series"]["times"])
        curve, top_edge, tri_edge = [], [], []
        for i in range(n_t):
            diffs, top_vals, tri_vals = [], [], []
            for g in pairs:
                top = g["chiral_topological"]["time_series"]
                tri = g["chiral_trivial"]["time_series"]
                top_vals.append(top["negativity_edge"][i])
                tri_vals.append(tri["negativity_edge"][i])
                diffs.append(top["negativity_edge"][i] - tri["negativity_edge"][i])
            curve.append(stats(diffs))
            top_edge.append(stats(top_vals))
            tri_edge.append(stats(tri_vals))
        means = [c["mean"] for c in curve]
        out[f"{profile}_d{disorder:g}"] = {
            "fraction_of_period": next(iter(pairs))["chiral_topological"]["time_series"]["fraction_of_period"],
            "contrast": curve,
            "negativity_edge_topological": top_edge,
            "negativity_edge_trivial": tri_edge,
            "mean_curve": means,
            "contrast_at_period_end": curve[-1],
            "sign_stable_over_period": bool(all(m > 0 for m in means) or all(m < 0 for m in means)),
            "minimum_mean": float(min(means)), "maximum_mean": float(max(means)),
        }
    return {"kind": "micromotion", "n_cells": int(rows[0]["n_cells"]),
            "n_records": len(rows), "conditions": out}


def summarize_chiralloss(path: Path) -> dict:
    """Paired sector contrasts and localization for the R5 loss-profile scan."""
    rows = load_records(path)
    groups = condition_groups(rows)
    out = {}
    for profile, disorder in sorted({(p, d) for (p, d, _) in groups}):
        pairs = [g for (p, d, _), g in groups.items()
                 if p == profile and d == disorder and set(g) == set(SECTORS)]
        if not pairs:
            continue
        contrast = stats(g["chiral_topological"]["negativity_edge"]
                         - g["chiral_trivial"]["negativity_edge"] for g in pairs)
        entry = {
            "paired_contrast": contrast,
            "negativity_edge_topological": stats(g["chiral_topological"]["negativity_edge"] for g in pairs),
            "negativity_edge_trivial": stats(g["chiral_trivial"]["negativity_edge"] for g in pairs),
            "edge_minus_center_topological": stats(
                g["chiral_topological"]["negativity_edge"] - g["chiral_topological"]["negativity_center"]
                for g in pairs),
        }
        if profile == "chiral_even_loss":
            sample_row = next(r for r in rows if r["loss_profile"] == profile)
            entry["odd_sublattice_rates_zero"] = bool(
                all(x == 0.0 for x in sample_row["rates"][1::2]))
        out[f"{profile}_d{disorder:g}"] = entry
    return {"kind": "chiralloss", "n_cells": int(rows[0]["n_cells"]),
            "n_records": len(rows), "conditions": out}


def sweep_observables(rows: list[dict], disorder: float, sample: int) -> dict:
    row = next(r for r in rows if float(r["disorder_level"]) == disorder
                and int(r["sample"]) == sample)
    return {
        "negativity_edge": row["negativity_edge"],
        "negativity_center": row["negativity_center"],
        "edge_minus_center": row["negativity_edge"] - row["negativity_center"],
    }


def summarize_sweep(path: Path) -> dict:
    """Per-intra statistics of the R4 sweep (one record set per intra value)."""
    rows = load_records(path)
    meta = {k: rows[0][k] for k in (
        "intra_A", "inter_B", "drive_ratio", "nu_0", "nu_pi",
        "winding_frame_A", "winding_frame_B", "resolution_stable",
        "bulk_point_candidate_pass", "minimum_zero_gap", "minimum_pi_gap")}
    per_disorder = {}
    for disorder in sorted({float(r["disorder_level"]) for r in rows}):
        block = {}
        for obs in ("negativity_edge", "negativity_center", "edge_minus_center"):
            block[obs] = stats(sweep_observables(rows, disorder, s)[obs]
                               for s in sorted({int(r["sample"]) for r in rows}))
        per_disorder[f"d{disorder:g}"] = block
    pooled = {obs: stats(sweep_observables(rows, d, s)[obs]
                         for d in sorted({float(r["disorder_level"]) for r in rows})
                         for s in sorted({int(r["sample"]) for r in rows}))
              for obs in ("negativity_edge", "negativity_center", "edge_minus_center")}
    return {"kind": "sweep", "n_cells": int(rows[0]["n_cells"]),
            "n_records": len(rows), "point": meta,
            "per_disorder": per_disorder, "pooled": pooled}


def paired_curvature(cells: int, center: float, left: float, right: float,
                     observable: str) -> dict:
    """Residual of the critical point against linear interpolation of its brackets.

    All three sweep points share the same disorder draws and trajectory seeds,
    so the residual is formed *within* a draw before aggregation.  A jump (or a
    kink) of the observable at the critical intra value shows up as a residual
    that differs from zero; a smooth curve gives a residual consistent with
    zero at the level of the pointwise errors.
    """

    def load(intra):
        path = SCAN_DIR / f"sweep_sweep_intra{intra:g}_{cells}cells.json"
        rows = load_records(path)
        return {(float(r["disorder_level"]), int(r["sample"])):
                sweep_observables(rows, float(r["disorder_level"]), int(r["sample"]))[observable]
                for r in rows}

    c, a, b = load(center), load(left), load(right)
    keys = sorted(set(c) & set(a) & set(b))
    if not keys:
        raise FileNotFoundError(f"no matched draws around intra={center}")
    #  x(c) - [x(a) + x(b)]/2  is the deviation from linear interpolation.
    residual = stats(c[k] - 0.5 * (a[k] + b[k]) for k in keys)
    return {
        "observable": observable, "center_intra": center,
        "left_intra": left, "right_intra": right, "n_pairs": len(keys),
        "center": stats(c[k] for k in keys),
        "left": stats(a[k] for k in keys),
        "right": stats(b[k] for k in keys),
        "linear_interpolation": stats(0.5 * (a[k] + b[k]) for k in keys),
        "residual_from_linear": residual,
        "z_score": (float(residual["mean"] / residual["se"])
                    if residual["se"] else 0.0),
    }


def main():
    out = {"status": "REVISION_SCANS_SUMMARIZED", "scans": {}}

    gamma_paths = (sorted(SCAN_DIR.glob("gamma_ratio_*_2cells.json"))
                   + sorted(SCAN_DIR.glob("gamma_matched_*_2cells.json")))
    for path in gamma_paths:
        matched = "_matched_" in path.name
        ratio = float(path.name.split("_")[2])
        rows = load_records(path)
        design = (json.loads(path.read_text())["design"]
                  if path.is_file() else {"gamma_edge": rows[0]["rates"][0],
                                          "gamma_bulk": rows[0]["rates"][1]})
        profiles = sorted({r["loss_profile"] for r in rows})
        out["scans"][f"gamma_{'matched' if matched else 'ratio'}_{ratio:g}"] = {
            "kind": "gamma_ratio_matched" if matched else "gamma_ratio",
            "gamma_edge": design["gamma_edge"],
            "gamma_bulk": design["gamma_bulk"],
            "budget_matched": bool(design.get("budget_matched", matched)),
            "budget_total": design.get("budget_total"),
            "n_cells": 2,
            "edge_low_loss": ({f"d{float(d):g}": summarize(rows, 2, "edge_low_loss", float(d))
                               for d in (0.0, 0.5)} if "edge_low_loss" in profiles else {}),
            "uniform_loss": ({f"d{float(d):g}": summarize(rows, 2, "uniform_loss", float(d))
                              for d in (0.0, 0.5)} if "uniform_loss" in profiles else {}),
        }

    for path in sorted(SCAN_DIR.glob("multiperiod_*T_*cells.json")):
        name = path.stem
        n_cells = int(name.split("_")[-1].replace("cells", ""))
        n_periods = int(name.split("_")[1].replace("T", ""))
        rows = load_records(path)
        out["scans"][name] = {
            "kind": "multiperiod",
            "n_periods": n_periods,
            "n_cells": n_cells,
            "edge_low_loss": {f"d{float(d):g}": summarize(rows, n_cells, "edge_low_loss", float(d)) for d in (0.0, 0.5)},
            "uniform_loss": {f"d{float(d):g}": summarize(rows, n_cells, "uniform_loss", float(d)) for d in (0.0, 0.5)},
        }

    nmax_paths = sorted(set(SCAN_DIR.glob("nmax_*_*cells.json"))
                        | set(SCAN_DIR.glob("nmax_*_*cells_nmc*.json")))
    for path in nmax_paths:
        name = path.stem
        cells_token = next(tok for tok in name.split("_") if tok.endswith("cells"))
        n_cells = int(cells_token.replace("cells", ""))
        n_max = int(name.split("_")[1])
        rows = load_records(path)
        out["scans"][name] = {
            "kind": "nmax",
            "n_max": n_max,
            "n_cells": n_cells,
            "sample_raised_extension": name.endswith("nmc200"),
            "edge_low_loss": {f"d{float(d):g}": summarize(rows, n_cells, "edge_low_loss", float(d)) for d in (0.0, 0.35)},
            "uniform_loss": {f"d{float(d):g}": summarize(rows, n_cells, "uniform_loss", float(d)) for d in (0.0, 0.35)},
        }

    for path in sorted(SCAN_DIR.glob("micromotion_*cells.json")):
        out["scans"][path.stem] = summarize_micromotion(path)

    for path in sorted(SCAN_DIR.glob("chiralloss_*cells.json")):
        out["scans"][path.stem] = summarize_chiralloss(path)

    sweep_paths = sorted(SCAN_DIR.glob("sweep_sweep_intra*_*cells.json"))
    for path in sweep_paths:
        summary = summarize_sweep(path)
        cells = summary["n_cells"]
        key = f"sweep_{cells}cells"
        out["scans"].setdefault(key, {"kind": "sweep", "n_cells": cells, "points": []})
        out["scans"][key].setdefault("points", []).append(summary)
    for key, block in out["scans"].items():
        if block.get("kind") == "sweep":
            block["points"].sort(key=lambda p: p["point"]["intra_A"])
            cells = block["n_cells"]
            present = [p["point"]["intra_A"] for p in block["points"]]
            tests = {}
            for label, center, left, right in (("pi_left", 0.15, 0.13, 0.17),
                                               ("zero", 0.65, 0.625, 0.675),
                                               ("pi_right", 1.15, 1.13, 1.17)):
                if center in present and left in present and right in present:
                    tests[label] = {
                        obs: paired_curvature(cells, center, left, right, obs)
                        for obs in ("negativity_edge", "edge_minus_center")}
            block["paired_smoothness_tests"] = tests

    out_path = ROOT / "data" / "revision_scans_summary.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
