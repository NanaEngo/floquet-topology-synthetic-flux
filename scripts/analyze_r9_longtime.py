#!/usr/bin/env python3
"""R9 analysis: long-time sign persistence and decay of the edge contrast.

Implements audit run R9 (AUDIT_RESEARCH_REFINEMENTS_20261001.md line 127).
For each (size, profile, disorder) stratum the paired topological-minus-trivial
edge negativity contrast is computed on the stroboscopic grid t_k = k T,
k = 0..64, and the analysis reports:

  * the contrast at the production snapshot t = T (sanity: must reproduce the
    production mean within the combined standard error of the overlapping
    draw subset, samples 0..199);
  * the quasi-stationary plateau value (median of t = 4T..16T) with its sign;
  * sign persistence: the set of k where the mean contrast keeps the sign it
    has at t = T, and the first k where the sign flips (if any);
  * the stationary time: the first k from which the contrast stays within
    2 standard errors of zero for the rest of the horizon;
  * an effective envelope decay rate (fit of log|mean contrast| over the
    decay window) as a size-resolved gap proxy.

Outputs:
    data/r9_longtime_summary.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

FILES = {2: ROOT / "data" / "r9_longtime" / "r9_longtime_2cells_NMC200.json",
         3: ROOT / "data" / "r9_longtime" / "r9_longtime_3cells_NMC200.json"}
PROD_SUMMARY = ROOT / "data" / "chiral_production_summary.csv"
OBS = "negativity_edge_topological_minus_trivial"
PLATEAU = (4, 16)  # stroboscopic indices treated as quasi-stationary plateau


def production_mean(cells: int, disorder: float, profile: str):
    """Production contrast mean and SE from the summary CSV."""
    import csv
    with PROD_SUMMARY.open() as fh:
        for row in csv.DictReader(fh):
            if (row["observable"] == OBS and int(row["n_cells"]) == cells
                    and float(row["disorder_level"]) == disorder
                    and row["loss_profile"] == profile):
                return float(row["mean"]), float(row["se"])
    return None, None


def paired_series(records: list, disorder: float, profile: str):
    """Paired contrast time series per sample.

    Args:
        records (list): R9 record rows.
        disorder (float): Disorder level filter (dimensionless).
        profile (str): Loss profile filter.

    Returns:
        ndarray: Array of shape (n_samples, n_times) of contrasts
            (dimensionless).
    """
    sel = [r for r in records if r["disorder_level"] == disorder
           and r["loss_profile"] == profile]
    topo = {r["sample"]: np.asarray(r["time_series"]["negativity_edge"])
            for r in sel if r["sector"] == "chiral_topological"}
    triv = {r["sample"]: np.asarray(r["time_series"]["negativity_edge"])
            for r in sel if r["sector"] == "chiral_trivial"}
    common = sorted(set(topo) & set(triv))
    return np.array([topo[s] - triv[s] for s in common])


def analyse_stratum(mat: np.ndarray) -> dict:
    """Mean/SE series, plateau, sign persistence, stationary time, decay rate.

    Args:
        mat (ndarray): (n_samples, n_times) paired contrasts (dimensionless).

    Returns:
        dict: Stratum diagnostics (dimensionless / per period).
    """
    mean = mat.mean(axis=0)
    se = mat.std(axis=0, ddof=1) / np.sqrt(mat.shape[0])
    k1 = 1  # t = T

    plateau_slice = slice(*PLATEAU)
    plateau = float(np.median(mean[plateau_slice]))
    plateau_se = float(np.median(se[plateau_slice]))

    sign_t = int(np.sign(mean[k1]))
    same_sign = [int(k) for k in range(1, len(mean))
                 if np.sign(mean[k]) == sign_t and abs(mean[k]) > 0]
    flip_k = next((k for k in range(k1 + 1, len(mean))
                   if np.sign(mean[k]) not in (0, sign_t)
                   and abs(mean[k]) > 2 * se[k]), None)

    stationary_k = None
    for k0 in range(1, len(mean)):
        if all(abs(mean[k]) <= 2 * se[k] for k in range(k0, len(mean))):
            stationary_k = k0
            break

    decay = {"rate_per_period": None, "tau_periods": None, "fit_window": None}
    if flip_k is None and stationary_k is not None and stationary_k > PLATEAU[1]:
        ks = np.arange(PLATEAU[1] + 1, stationary_k)
        vals = np.abs(mean[ks])
        ok = vals > 0
        if np.sum(ok) >= 4:
            slope = float(np.polyfit(ks[ok], np.log(vals[ok]), 1)[0])
            if slope < 0:
                decay = {"rate_per_period": -slope,
                         "tau_periods": -1.0 / slope,
                         "fit_window": [int(ks[ok][0]), int(ks[ok][-1])]}

    return {
        "n_samples": int(mat.shape[0]),
        "contrast_at_T": {"mean": float(mean[k1]), "se": float(se[k1]),
                          "ci95": [float(mean[k1] - 1.96 * se[k1]),
                                   float(mean[k1] + 1.96 * se[k1])]},
        "plateau_median_4T_16T": {"mean": plateau, "se": plateau_se,
                                  "sign": int(np.sign(plateau))},
        "sign_at_T": sign_t,
        "sign_persistence_first_last": [same_sign[0], same_sign[-1]] if same_sign else [],
        "sign_flip_first_k": flip_k,
        "stationary_first_k": stationary_k,
        "contrast_at_horizon": {"mean": float(mean[-1]), "se": float(se[-1])},
        "envelope_decay": decay,
        "mean_series": [float(v) for v in mean],
        "se_series": [float(v) for v in se],
    }


def main() -> None:
    """Aggregate R9 series for all strata and write the summary JSON."""
    summary = {}
    for cells, path in FILES.items():
        data = json.loads(path.read_text())
        records = data["records"]
        for disorder in sorted({r["disorder_level"] for r in records}):
            for profile in sorted({r["loss_profile"] for r in records}):
                mat = paired_series(records, disorder, profile)
                entry = analyse_stratum(mat)
                pm, pse = production_mean(cells, disorder, profile)
                if pm is not None:
                    # production mean covers samples 0..499; R9 covers 0..199.
                    # consistency check within the combined SE (loose but
                    # fail-closed if the sign or scale disagrees).
                    combined = float(np.hypot(entry["contrast_at_T"]["se"], pse))
                    entry["production_cross_check"] = {
                        "production_mean": pm, "production_se": pse,
                        "delta": entry["contrast_at_T"]["mean"] - pm,
                        "within_3_combined_se": bool(
                            abs(entry["contrast_at_T"]["mean"] - pm) <= 3 * combined),
                    }
                summary[f"{cells}cells_d{disorder:g}_{profile}"] = entry

    out = {
        "status": "COMPLETED_R9_LONGTIME_ANALYSIS",
        "audit_source": "AUDIT_RESEARCH_REFINEMENTS_20261001.md (line 127)",
        "time_grid": "stroboscopic t_k = k*PERIOD, k = 0..64",
        "plateau_window": list(PLATEAU),
        "strata": summary,
    }
    dest = ROOT / "data" / "r9_longtime_summary.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dest}")
    for key, e in summary.items():
        t = e["contrast_at_T"]
        p = e["plateau_median_4T_16T"]
        pc = e.get("production_cross_check", {})
        print(f"  {key}: d(T)={t['mean']:+.5f}±{t['se']:.5f} "
              f"plateau={p['mean']:+.5f} sign_flip_k={e['sign_flip_first_k']} "
              f"stationary_k={e['stationary_first_k']} "
              f"tau={e['envelope_decay']['tau_periods']} "
              f"prod_xcheck={pc.get('within_3_combined_se')}")


if __name__ == "__main__":
    main()
