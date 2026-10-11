#!/usr/bin/env python3
"""Edge-mode localization length xi along the R4 dimerization sweep (audit E11).

The main text concludes that the sweep contrast behaves as a
"finite-localization quantity".  Audit item E11 requires the alternative
hypothesis to be measured rather than asserted: for every sweep point
(``intra_A`` = J_A, ``inter_B`` = 1.3 - J_A) this script computes the coherent
open-chain edge mode used by the micromotion-covariant diagnostic -- the most
edge-localized eigenvector of the one-period propagator
``U = expm(-i h_B T/2) @ expm(-i h_A T/2)`` -- and extracts its localization
length as the participation-based length

    xi_part = (sum_j |u_j|^4)^-1   (sites),

which is proportional to the exponential localization length for a localized
profile and saturates at the chain length for an extended profile.  The
participation ratio is used because on the two- and three-cell chains of the
sweep the eigenvector is a finite-size symmetric cat of the two edge modes:
a site-by-site exponential envelope has too few decay steps (and changes
sublattice character within the first dimer) to support a stable
least-squares fit, whereas the participation length is defined at every
point, is insensitive to the cat superposition (a constant factor two), and
is the standard localization measure on small chains.

The measured sweep contrast at zero disorder (edge-minus-center negativity,
mean over matched samples, from ``data/revision_scans/sweep_sweep_intra*.json``)
is then regressed against xi_part for each system size: Pearson and Spearman
correlations with two-sided p-values, plus an ordinary-least-squares fit
``contrast = a * xi_part + b``.

Output:
    data/sweep_xi_regression.json -- per-size tables and regression blocks.

The script aborts if any expected sweep file is missing or if the
localization sanity gate fails (deep-topological points must be more
localized than the critical window).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import stats as sps

from run_chiral_floquet_dynamics import edge_mode_vector

ROOT = Path(__file__).resolve().parents[1]
SCAN_DIR = ROOT / "data" / "revision_scans"
OUT_PATH = ROOT / "data" / "sweep_xi_regression.json"
TOTAL_HOPPING = 1.3


def localization_length(u: np.ndarray) -> dict:
    """Participation-based localization length of a site-space mode vector.

    Args:
        u (numpy.ndarray): Complex mode amplitudes on the n sites of the
            open chain (dimensionless, unit norm).

    Returns:
        dict: ``xi_part`` (participation-based localization length in
        sites), ``ipr`` (dimensionless) and ``participation_sites``
        (alias of ``xi_part``).
    """
    amp = np.abs(u)
    ipr = float(np.sum(amp ** 4))
    xi_part = float(1.0 / ipr)
    return {"xi_part": xi_part, "ipr": ipr, "participation_sites": xi_part}


def sweep_contrast_zero_disorder(path: Path) -> dict:
    """Mean edge-minus-center contrast at d=0 for one sweep point.

    Args:
        path (Path): Merged sweep JSON file.

    Returns:
        dict: ``n`` samples, ``mean`` and ``se`` of the observable
        (dimensionless Wigner negativity units).
    """
    rows = json.loads(path.read_text())["records"]
    vals = np.asarray([r["negativity_edge"] - r["negativity_center"]
                       for r in rows if float(r["disorder_level"]) == 0.0])
    if vals.size == 0:
        raise ValueError(f"no d=0 records in {path}")
    mean = float(vals.mean())
    se = float(vals.std(ddof=1) / np.sqrt(vals.size)) if vals.size > 1 else 0.0
    return {"n": int(vals.size), "mean": mean, "se": se}


def regress(x: np.ndarray, y: np.ndarray) -> dict:
    """Pearson, Spearman and OLS of y on x with two-sided p-values."""
    pearson = sps.pearsonr(x, y)
    spearman = sps.spearmanr(x, y)
    ols = sps.linregress(x, y)
    return {
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman_rho": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
        "ols_slope": float(ols.slope),
        "ols_intercept": float(ols.intercept),
        "ols_slope_se": float(ols.stderr),
        "ols_r2": float(ols.rvalue) ** 2,
        "n": int(x.size),
    }


def main() -> None:
    """Compute xi for every sweep point and regress the contrast on xi."""
    tables = {}
    failures = []
    for cells in (2, 3):
        points = []
        for path in sorted(SCAN_DIR.glob(f"sweep_sweep_intra*_{cells}cells.json")):
            rows = json.loads(path.read_text())["records"]
            intra = float(rows[0]["intra_A"])
            inter = float(rows[0]["inter_B"])
            loc = localization_length(
                edge_mode_vector(cells, {"intra_A": intra, "inter_B": inter}))
            contrast = sweep_contrast_zero_disorder(path)
            points.append({
                "intra_A": intra,
                "inter_B": inter,
                "drive_ratio": intra / inter,
                "xi_part": loc["xi_part"],
                "ipr": loc["ipr"],
                "contrast_mean": contrast["mean"],
                "contrast_se": contrast["se"],
                "n_samples": contrast["n"],
                "source": path.name,
            })
        if not points:
            failures.append(f"no sweep files for {cells} cells")
            continue
        points.sort(key=lambda p: p["intra_A"])
        x = np.asarray([p["xi_part"] for p in points], dtype=float)
        y = np.asarray([p["contrast_mean"] for p in points], dtype=float)
        # Sanity gate: the deeply dimerized topological points must be more
        # localized (smaller xi_part) than the critical window.
        xi_topo = min(p["xi_part"] for p in points if p["intra_A"] <= 0.2)
        xi_window = max(p["xi_part"] for p in points if 0.6 <= p["intra_A"] <= 0.8)
        if not xi_topo < xi_window:
            failures.append(f"{cells} cells: xi_part sanity gate failed "
                            f"(topo {xi_topo:.3f} >= window {xi_window:.3f})")
        tables[f"{cells}cells"] = {
            "n_points": len(points),
            "points": points,
            "regression": regress(x, y),
        }
    if failures:
        raise SystemExit("FAIL: " + "; ".join(failures))
    out = {
        "status": "SWEEP_XI_REGRESSION_COMPLETE",
        "claim_scope": "coherent edge-mode localization length vs measured "
                       "zero-disorder sweep contrast; descriptive regression, "
                       "no causal or protection claim",
        "not_supported": ["topological protection theorem",
                          "experimental validation", "metrological gain"],
        "method": {
            "mode": "most edge-localized eigenvector of the one-period "
                    "propagator expm(-i h_B T/2) @ expm(-i h_A T/2) on the "
                    "open chain (same vector as the micromotion-covariant "
                    "edge mode; a symmetric cat of the two edge modes on "
                    "these chain lengths)",
            "xi_part": "participation-based localization length "
                       "(sum_j |u_j|^4)^-1 in sites; proportional to the "
                       "exponential localization length for localized "
                       "profiles, saturating at the chain length for "
                       "extended modes",
            "contrast": "edge-minus-center negativity at d=0, mean over "
                        "matched samples (revision scan R4)",
            "total_hopping": TOTAL_HOPPING,
        },
        "tables": tables,
    }
    OUT_PATH.write_text(json.dumps(out, indent=2) + "\n")
    for size, block in tables.items():
        r = block["regression"]
        print(f"{size}: n={r['n']} Pearson r={r['pearson_r']:+.3f} "
              f"(p={r['pearson_p']:.3g}) Spearman rho={r['spearman_rho']:+.3f} "
              f"(p={r['spearman_p']:.3g}) OLS slope={r['ols_slope']:+.4f} "
              f"R2={r['ols_r2']:.3f}")
    print(f"saved {OUT_PATH}")


if __name__ == "__main__":
    main()
