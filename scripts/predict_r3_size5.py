"""Frozen out-of-sample prediction of the 5-cell (10-mode) edge contrast.

Pre-registers, BEFORE the R3 size-5 production run is submitted, a prediction
of the edge topological-minus-trivial negativity contrast at five cells from
the three existing production sizes (2, 3, 4 cells). This implements the L3
out-of-sample requirement of ``AUDIT_RESEARCH_REFINEMENTS_20261001.md``:

    "fit contrast with diagnostics (edge participation, Gamma_e exposure,
     gamma ratio), predict 5th size out-of-sample"

Two pre-registered models are fitted to the three production means per
(disorder, loss-profile) cell:

    M1: Delta(N) = Delta_inf + a * exp(-(N-2)/xi)          (3 parameters)
    M2: Delta(N) = Delta_inf + a / N                        (2 parameters)

The interval combines the model spread with the production standard errors
at 5 cells (rescaled by the ratio of sample sizes if needed). The resulting
JSON is written once and must not be rewritten after the R3 run starts; the
post-run analysis scores the measured value against this file.

The size-dependent diagnostics are the single-particle edge-mode inverse
participation ratio (IPR) and the loss-rate exposure
Gamma_e = sum_j gamma_j |v_j|^2 / sum_j |v_j|^2 of the edge mode under the
edge-low-loss profile (edge sites 1e-4, bulk 0.2, dimensionless rates in
units of the drive frequency), together with the constant bulk/edge rate
ratio.

Outputs:
    data/r3_size5_prediction.json

Run (before ``sbatch`` submission of the R3 run):
    mamba run -n qom python scripts/predict_r3_size5.py
"""

import csv
import json
import time
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_chiral_floquet_dynamics import CANDIDATES, edge_mode_vector  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "data" / "chiral_production_summary.csv"
OUTPUT = ROOT / "data" / "r3_size5_prediction.json"

OBSERVABLE = "negativity_edge_topological_minus_trivial"
SIZES = (2, 3, 4)
TARGET_SIZE = 5
EXPECTED_R3_N = 300  # R3 audit spec: N_MC >= 300
GAMMA_EDGE = 1e-4  # dimensionless, units of drive frequency
GAMMA_BULK = 0.2   # dimensionless, units of drive frequency


def load_contrasts():
    """Read production contrast means and standard errors per condition.

    Returns:
        dict: {(disorder, profile, n_cells): {"mean": float, "se": float, "n": int}}
            with means of the edge topological-minus-trivial negativity
            contrast (dimensionless) from the production summary CSV.
    """
    out = {}
    with SUMMARY.open() as fh:
        for row in csv.DictReader(fh):
            if row["observable"] != OBSERVABLE:
                continue
            key = (float(row["disorder_level"]), row["loss_profile"], int(row["n_cells"]))
            out[key] = {"mean": float(row["mean"]), "se": float(row["se"]),
                        "n": int(row["n"]), "ci95_low": float(row["ci95_low"]),
                        "ci95_high": float(row["ci95_high"])}
    return out


def fit_models(sizes, means):
    """Fit the two pre-registered saturation models to three size points.

    Args:
        sizes (array-like): Cell counts (dimensionless), length 3.
        means (array-like): Mean contrasts (dimensionless), length 3.

    Returns:
        dict: Model predictions at ``TARGET_SIZE`` with coefficients.
    """
    sizes = np.asarray(sizes, float)
    means = np.asarray(means, float)

    # M1: three-parameter exponential saturation, exact fit of three points.
    # Solve linearly in (Delta_inf, a) for each decay factor b = exp(-1/xi):
    #   means[k] = D + a * b**(sizes[k]-2)
    best = None
    for b in b_grid(means):
        Amat = np.column_stack([np.ones(3), b ** (sizes - 2)])
        coef, *_ = np.linalg.lstsq(Amat, means, rcond=None)
        resid = float(np.max(np.abs(Amat @ coef - means)))
        if best is None or resid < best[0]:
            best = (resid, b, coef)
    resid, b, (D_inf, a) = best
    xi = -1.0 / np.log(b) if 0 < b < 1 else float("inf")
    pred_m1 = float(D_inf + a * b ** (TARGET_SIZE - 2))

    # M2: two-parameter 1/N saturation, least squares on three points.
    A2 = np.column_stack([np.ones(3), 1.0 / sizes])
    (D2, a2), *_ = np.linalg.lstsq(A2, means, rcond=None)
    pred_m2 = float(D2 + a2 / TARGET_SIZE)

    return {
        "M1_exponential": {"pred": pred_m1, "Delta_inf": float(D_inf),
                           "a": float(a), "xi_cells": float(xi),
                           "max_residual": resid},
        "M2_inverse_N": {"pred": pred_m2, "Delta_inf": float(D2), "a": float(a2)},
    }


def b_grid(means):
    """Candidate decay factors b for the exponential model.

    Args:
        means (array-like): Three contrast means (dimensionless).

    Returns:
        ndarray: Log-spaced grid of b in (0.001, 0.999).
    """
    return np.linspace(0.001, 0.999, 999)


def edge_diagnostics(n_cells):
    """Compute size-dependent edge-mode diagnostics.

    Args:
        n_cells (int): Number of unit cells (dimensionless).

    Returns:
        dict: IPR (dimensionless), Gamma_e exposure (units of drive
            frequency), and the constant bulk/edge rate ratio.
    """
    v = np.abs(np.asarray(edge_mode_vector(n_cells, CANDIDATES["chiral_topological"]), complex))
    w = v ** 2
    w = w / w.sum()
    ipr = float(np.sum(w ** 2))
    # site layout: within cell (A, B); edges are site 0 and site n_modes-1
    n_modes = len(v)
    gamma_vec = np.full(n_modes, GAMMA_BULK)
    gamma_vec[0] = GAMMA_EDGE
    gamma_vec[-1] = GAMMA_EDGE
    gamma_e = float(np.sum(gamma_vec * w))
    return {"ipr": ipr, "gamma_e_exposure": gamma_e,
            "gamma_bulk_over_edge": GAMMA_BULK / GAMMA_EDGE}


def main():
    """Fit, freeze and write the out-of-sample prediction JSON."""
    contrasts = load_contrasts()
    conditions = sorted({(d, p) for (d, p, n) in contrasts if n in SIZES})

    predictions = {}
    for disorder, profile in conditions:
        try:
            means = [contrasts[(disorder, profile, n)]["mean"] for n in SIZES]
            ses = [contrasts[(disorder, profile, n)]["se"] for n in SIZES]
            ns = [contrasts[(disorder, profile, n)]["n"] for n in SIZES]
        except KeyError:
            continue
        models = fit_models(SIZES, means)
        preds = [m["pred"] for m in models.values()]
        # measurement SE of the future R3 mean: pooled production sd / sqrt(N_R3)
        sd_pooled = float(np.mean([se * np.sqrt(n) for se, n in zip(ses, ns)]))
        se_target = sd_pooled / np.sqrt(EXPECTED_R3_N)
        # pre-registered interval: model envelope (all models must agree on the
        # sign) widened by the measurement SE of the future 5-cell mean
        lo = min(preds) - 1.96 * se_target
        hi = max(preds) + 1.96 * se_target
        center = float(np.mean(preds))
        predictions[f"disorder_{disorder}_{profile}"] = {
            "disorder_level": disorder, "loss_profile": profile,
            "observed_means": dict(zip(map(str, SIZES), means)),
            "models": models,
            "ensemble_prediction": {"mean": center, "ci95_low": lo, "ci95_high": hi,
                                    "se_target": se_target},
            "sign_prediction": "positive" if lo > 0 else (
                "negative" if hi < 0 else "unresolved"),
        }

    diagnostics = {str(n): edge_diagnostics(n) for n in (*SIZES, TARGET_SIZE)}

    out = {
        "status": "FROZEN_BEFORE_R3_RUN",
        "purpose": "L3 out-of-sample prediction of the 5-cell edge contrast",
        "source": "data/chiral_production_summary.csv (N_MC=500 production)",
        "audit_ref": "AUDIT_RESEARCH_REFINEMENTS_20261001.md (line 155, L3)",
        "target_size_cells": TARGET_SIZE,
        "target_n_modes": 2 * TARGET_SIZE,
        "pre_registered_models": ["M1_exponential", "M2_inverse_N"],
        "predictions": predictions,
        "edge_mode_diagnostics": diagnostics,
        "created_unix": time.time(),
        "created_iso": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "scoring_rule": ("data/r3_size5_analysis.json must score the measured "
                         "5-cell mean against ensemble_prediction and "
                         "sign_prediction without rewriting this file"),
    }
    if OUTPUT.exists():
        raise SystemExit(f"refusing to overwrite frozen prediction {OUTPUT}")
    OUTPUT.write_text(json.dumps(out, indent=2) + "\n")
    print(f"froze {OUTPUT}")
    for k, v in predictions.items():
        e = v["ensemble_prediction"]
        print(f"  {k}: {e['mean']:+.6f} "
              f"[{e['ci95_low']:+.6f},{e['ci95_high']:+.6f}] "
              f"-> {v['sign_prediction']}")


if __name__ == "__main__":
    main()
