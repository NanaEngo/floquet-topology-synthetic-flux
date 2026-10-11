#!/usr/bin/env python3
"""Recompute the bulk-participation decay-bound verification (SI Sec. S9).

Deterministic producer for the two SI tables derived from the chiral production
files:

  * ``tab:bound-verification`` -- zero-disorder edge-low-loss, topological
    sector, mean final negativity:

        Gamma_meas = -ln[ delta_-(T) / delta_-(0) ] / T ,
        delta_-(0) = 0.19623 (cat-state edge-site initial value on the
        production 31x31 grid),

    with Gamma_eff = sum_j gamma_j |u_{e,j}|^2 using the record's own
    disorder-drawn rates and the clean topological edge-mode weights.

  * ``tab:extended-bound`` -- per-condition mean-level ratio (valid-only
    condition-level ratio, see ``table_extended_bound``) and per-realization
    compliance fraction Gamma_meas >= Gamma_eff, over the full grid
    3 sizes x 2 profiles x 4 disorder levels, topological-sector records only.

For every size the script uses the largest available N_MC file (same selection
rule as ``analyze_chiral_production._dedupe_by_cells``), so the tables follow
the current primary production automatically.

Output: data/participation_bound_verification.json

Usage:
    mamba run -n qom python scripts/verify_participation_bound.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from run_chiral_floquet_dynamics import CANDIDATES, edge_mode_vector

ROOT = Path(__file__).resolve().parents[1]
PROD = ROOT / "data" / "chiral_production"
OUT = ROOT / "data" / "participation_bound_verification.json"
PERIOD = 2.0 * np.pi
DELTA0 = 0.19623  # deterministic cat initial negativity on the 31x31 grid
PROFILES = ("edge_low_loss", "uniform_loss")
DISORDER_LEVELS = (0.0, 0.5, 1.0, 1.5)


def best_file(n_cells: int) -> Path:
    """Return the production file with the largest N_MC for ``n_cells``."""
    cands = list(PROD.glob(f"chiral_floquet_{n_cells}cells_NMC*.json"))
    if not cands:
        raise SystemExit(f"no production file for {n_cells} cells")
    return max(cands, key=lambda p: int(json.loads(p.read_text())["design"]["n_mc"]))


def edge_weights(n_cells: int) -> np.ndarray:
    """Clean topological-sector edge-mode weights |u_{e,j}|^2 for the size."""
    u = edge_mode_vector(n_cells, CANDIDATES["chiral_topological"])
    w = np.abs(u) ** 2
    return w / w.sum()


def load_records(path: Path):
    return json.loads(path.read_text())["records"]


def gamma_eff(rates: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(rates * weights))


def gamma_meas(neg_t: float) -> float | None:
    """-ln[delta_-(T)/delta_-(0)]/T per realization.

    A realization whose negativity has fully decayed (``neg_t <= 0``)
    has an unbounded measured decay rate, so the bound
    ``Gamma_meas >= Gamma_eff`` is trivially satisfied; we return
    ``+inf`` for those (counted as compliant, excluded from the
    mean-level ratio, which uses the condition-mean negativity).
    ``None`` is returned only for non-finite inputs.
    """
    if not np.isfinite(neg_t):
        return None
    if neg_t <= 0.0:
        return math.inf
    return -math.log(neg_t / DELTA0) / PERIOD


def table_bound_verification() -> list[dict]:
    rows = []
    for n_cells in (2, 3, 4):
        path = best_file(n_cells)
        weights = edge_weights(n_cells)
        recs = [
            r for r in load_records(path)
            if r["loss_profile"] == "edge_low_loss"
            and float(r["disorder_level"]) == 0.0
            and r["sector"] == "chiral_topological"
        ]
        negs = np.array([float(r["negativity_edge"]) for r in recs])
        g_eff = gamma_eff(np.asarray(
            [float(x) for x in recs[0]["rates"]], dtype=float), weights)
        g_meas = gamma_meas(float(negs.mean()))
        rows.append({
            "n_cells": n_cells,
            "n_modes": 2 * n_cells,
            "n_records": len(recs),
            "file": str(path.relative_to(ROOT)),
            "n_mc": int(json.loads(path.read_text())["design"]["n_mc"]),
            "edge_weight_site0": float(weights[0]),
            "gamma_eff": g_eff,
            "negativity_T_mean": float(negs.mean()),
            "negativity_0": DELTA0,
            "gamma_meas": float(g_meas),
            "ratio": float(g_meas / g_eff),
        })
    return rows


def table_extended_bound() -> dict:
    """Per-condition extended verification over the full 24-condition grid.

    Two statistics per (profile, disorder, size) condition:

    * ``mean_level_ratio`` -- valid-only condition-level ratio: Gamma_meas
      evaluated at the mean final negativity over the realizations that
      still carry a positive negativity at T, divided by the mean Gamma_eff
      over those same realizations;
      realizations whose negativity has fully decayed (delta_-(T) <= 0)
      have no finite Gamma_meas and are excluded from this statistic only.
      A per-realization average is deliberately not used: it is dominated by
      realizations whose Gamma_eff is driven toward zero by the rate clamping.
      Note that ``tab:bound-verification`` evaluates Gamma_meas at the mean
      over *all* realizations (its delta_-(T) column); the extended grid
      restricts the ratio to delta_-(T) > 0 realizations, which is why its
      upper range sits below the three-cell entry of that table.
    * ``compliance`` -- fraction of *all* realizations with
      Gamma_meas >= Gamma_eff, evaluated by the ordinary rule: a fully
      decayed realization (Gamma_meas = +inf) and a zero-rate realization
      (Gamma_eff = 0, bound reduces to delta_-(T) <= delta_-(0)) are both
      tested against that rule, never skipped. Only non-finite negativities
      leave the population.
    """
    out = {}
    for profile in PROFILES:
        for disorder in DISORDER_LEVELS:
            per_size = {}
            for n_cells in (2, 3, 4):
                path = best_file(n_cells)
                weights = edge_weights(n_cells)
                all_neg, valid_neg, valid_eff, comply = [], [], [], 0
                total = 0
                for r in load_records(path):
                    if r["loss_profile"] != profile:
                        continue
                    if float(r["disorder_level"]) != disorder:
                        continue
                    if r["sector"] != "chiral_topological":
                        continue
                    dt = float(r["negativity_edge"])
                    g_m = gamma_meas(dt)
                    if g_m is None:
                        continue
                    total += 1
                    all_neg.append(dt)
                    g_e = gamma_eff(
                        np.asarray([float(x) for x in r["rates"]], dtype=float),
                        weights,
                    )
                    if g_m >= g_e:
                        comply += 1
                    if dt > 0.0 and g_e > 0.0:
                        valid_neg.append(dt)
                        valid_eff.append(g_e)
                mean_neg = float(np.mean(valid_neg)) if valid_neg else 0.0
                mean_g_eff = float(np.mean(valid_eff)) if valid_eff else 0.0
                g_m_mean = gamma_meas(mean_neg)
                ratio = (g_m_mean / mean_g_eff
                         if g_m_mean is not None and mean_g_eff > 0
                         else float("nan"))
                per_size[n_cells] = {
                    "n_records": total,
                    "n_fully_decayed": sum(1 for t in all_neg if t <= 0.0),
                    "n_ratio_valid": len(valid_neg),
                    "mean_negativity_T": mean_neg,
                    "mean_level_ratio": float(ratio),
                    "compliance": comply / total if total else float("nan"),
                }
            rs = [v["mean_level_ratio"] for v in per_size.values()]
            cs = {k: v["compliance"] for k, v in per_size.items()}
            out[f"{profile}/d={disorder:g}"] = {
                "per_size": per_size,
                "ratio_range": [min(rs), max(rs)],
                "compliance": cs,
            }
    return out


def main() -> None:
    verification = table_bound_verification()
    extended = table_extended_bound()
    report = {
        "status": "PARTICIPATION_BOUND_RECOMPUTED",
        "scope": "empirical finite-time diagnostic; not an independently proved bound",
        "definition": {
            "delta_0": DELTA0,
            "gamma_meas": "-ln[delta_-(T)/delta_-(0)]/T over T=2*pi; a realization with delta_-(T) <= 0 has an unbounded measured rate (Gamma_meas = +inf)",
            "gamma_eff": "sum_j gamma_j |u_{e,j}|^2 with the record's own disorder-drawn rates and the clean topological edge-mode weights; rates are clamped at zero, so gamma_eff = 0 is possible at large disorder",
            "mean_level_ratio": "valid-only condition-level ratio: Gamma_meas at the mean negativity over realizations with delta_-(T) > 0, divided by the mean Gamma_eff over those same realizations (tab:bound-verification instead uses the mean over all realizations)",
            "compliance": "fraction of all finite-negativity realizations with Gamma_meas >= Gamma_eff, evaluated by the ordinary rule (fully decayed: +inf >= gamma_eff; zero-rate: Gamma_meas >= 0)",
            "file_selection": "largest available N_MC per size (same rule as analyze_chiral_production)",
            "sector": "chiral_topological records only",
        },
        "bound_verification": verification,
        "extended_bound": extended,
    }
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    for row in verification:
        print(f"{row['n_cells']} cells (N_MC={row['n_mc']}): "
              f"Gamma_meas={row['gamma_meas']:.4f} Gamma_eff={row['gamma_eff']:.4f} "
              f"ratio={row['ratio']:.2f}")
    for key, val in extended.items():
        print(f"{key}: ratio_range={val['ratio_range'][0]:.2f}-{val['ratio_range'][1]:.2f} "
              f"compliance={ {k: round(v * 100, 1) for k, v in val['compliance'].items()} }")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
