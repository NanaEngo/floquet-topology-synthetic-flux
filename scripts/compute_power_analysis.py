#!/usr/bin/env python3
"""Eight-mode statistical power statement (2026-08-17 adversarial audit, item 1.1).

Deterministic computation from the production 4-cell/8-mode data.  The
zero-disorder edge-low-loss paired contrasts (n = N_MC) are the independent
matched units.  We report the minimal |effect| a two-sided paired t-test can
resolve at alpha = 0.05, power = 0.8, for the attained sample size (100) and
for the size used at 2 and 3 cells (500):

    delta_min = (t_{1-alpha/2, n-1} + t_{1-beta, n-1}) * s / sqrt(n)

No distributional guarantee is asserted beyond the data at hand; this is a
descriptive power statement to bound the detectable effect size, not a
p-value test.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
ALPHA = 0.05
POWER = 0.80


def paired_contrasts(path: Path, profile: str, disorder: float) -> list[float]:
    data = json.loads(path.read_text())
    groups = {}
    for r in data["records"]:
        key = (r["loss_profile"], float(r["disorder_level"]), int(r["sample"]))
        groups.setdefault(key, {})[r["sector"]] = r
    out = []
    for (p, d, _), pair in groups.items():
        if p == profile and d == disorder:
            out.append(pair["chiral_topological"]["negativity_edge"]
                       - pair["chiral_trivial"]["negativity_edge"])
    return out


def detectable_effect(values: np.ndarray) -> float:
    n = len(values)
    s = float(np.std(values, ddof=1))
    t1 = stats.t.ppf(1 - ALPHA / 2, n - 1)
    t2 = stats.t.ppf(POWER, n - 1)
    return (t1 + t2) * s / math.sqrt(n)


def main():
    prod = ROOT / "data" / "chiral_production"
    files = sorted(prod.glob("chiral_floquet_*cells_NMC*.json"))
    if not files:
        raise SystemExit("no production files found")
    best: dict[int, tuple[int, Path]] = {}
    for f in files:
        design = json.loads(f.read_text())["design"]
        cells = int(design["n_cells"])
        n_mc = int(design["n_mc"])
        if cells not in best or n_mc > best[cells][0]:
            best[cells] = (n_mc, f)
    sizes = []
    for _, f in sorted(best.values(),
                       key=lambda t: int(json.loads(t[1].read_text())["design"]["n_cells"])):
        data = json.loads(f.read_text())
        cells = int(data["design"]["n_cells"])
        n_mc = int(data["design"]["n_mc"])
        contrasts = paired_contrasts(f, "edge_low_loss", 0.0)
        contrasts = np.asarray(contrasts)
        assert len(contrasts) == n_mc
        sizes.append({
            "cells": cells, "modes": 2 * cells, "n_mc": n_mc,
            "mean": float(contrasts.mean()),
            "sd": float(contrasts.std(ddof=1)),
            "detectable_effect_at_attained_n": detectable_effect(contrasts),
            "detectable_effect_at_n500": float(detectable_effect(contrasts) * math.sqrt(n_mc / 500)),
        })
    summary = {
        "status": "POWER_STATEMENT_COMPUTED",
        "scope": "descriptive power statement over the independent disorder units; no p-value testing",
        "design": {
            "test": "two-sided paired t-test",
            "alpha": ALPHA, "power": POWER,
            "observable": "zero-disorder edge-low-loss matched contrast (topological - trivial)",
            "independent_unit": "disorder realization",
            "formula": "delta_min = (t_{1-alpha/2,n-1} + t_{1-beta,n-1}) s / sqrt(n)",
        },
        "sizes": sizes,
        "software": {"numpy": np.__version__, "scipy": __import__("scipy").__version__},
    }
    out = ROOT / "data" / "power_analysis.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["sizes"], indent=2))


if __name__ == "__main__":
    main()
