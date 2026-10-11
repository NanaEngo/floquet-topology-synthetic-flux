#!/usr/bin/env python3
"""Trajectory-count extrapolation of the negativity contrast (audit A6 ladder).

The Monte Carlo production estimator evaluates the nonlinear Wigner-void
negativity on the trajectory-averaged state with ``ntraj = 10`` trajectories
per sample.  On conditions where the exact negativity is near zero this
estimator carries a positive finite-trajectory bias that differs between the
two sectors and therefore shifts the paired contrast
(``run_solver_contrast_check.py`` measures the bias directly at two and
three cells, where the density-matrix master equation is tractable).

This script assembles, per chain size, the contrast measured at trajectory
counts ``N = 10`` (archived production records, ``ntraj = 10``) and the
ladder points ``N`` produced by ``run_solver_contrast_check.py`` at
``ntraj > 10``, and fits the bias model

    mean(N) = C + c/N + d/N^2          (model B, unweighted, 3 parameters)
    mean(N) = C + c/N                  (model A, inverse-variance weighted)

where ``C`` is the bias-free (infinite-trajectory) contrast -- the exact
master-equation value.  At two and three cells the exact ``mesolve`` contrast
is available, so each model is validated against it before it is applied to
four and five cells, where the density-matrix solve is infeasible.

Output:
    data/solver_contrast_ladder.json

The script aborts if a required input is missing or if the preferred model
fails validation at both exactly solvable sizes.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CHECK_DIR = ROOT / "data"
PROD_DIR = ROOT / "data" / "chiral_production"
OUT_PATH = ROOT / "data" / "solver_contrast_ladder.json"
PROFILE = "edge_low_loss"
DISORDER = 0.0
VALIDATION_TOL = 0.0005  # absolute contrast units


def archive_contrast(n_cells: int) -> dict:
    """Production N=10 contrast statistics at zero disorder, edge profile.

    Args:
        n_cells (int): Chain size (unit cells).

    Returns:
        dict: ``ntraj``, ``n`` samples, ``mean``/``sd``/``se`` of the paired
        topological-minus-trivial negativity contrast (dimensionless) and the
        production file name.
    """
    paths = sorted(PROD_DIR.glob(f"chiral_floquet_{n_cells}cells_NMC*.json"),
                   key=lambda p: int(p.stem.rsplit("NMC", 1)[1]),
                   reverse=True)
    if not paths:
        raise FileNotFoundError(f"no production file for {n_cells} cells")
    path = paths[0]
    rows = json.loads(path.read_text())["records"]
    by_key = {(r["sample"], r["sector"]): r for r in rows
              if r["loss_profile"] == PROFILE
              and float(r["disorder_level"]) == DISORDER}
    samples = sorted({s for s, sec in by_key if sec == "chiral_topological"})
    vals = np.asarray([by_key[(s, "chiral_topological")]["negativity_edge"]
                       - by_key[(s, "chiral_trivial")]["negativity_edge"]
                       for s in samples])
    sd = float(vals.std(ddof=1))
    return {"ntraj": 10, "n": int(vals.size), "mean": float(vals.mean()),
            "sd": sd, "se": float(sd / np.sqrt(vals.size)),
            "source": path.name}


def ladder_points(n_cells: int) -> tuple[list[dict], dict | None]:
    """Ladder contrast points and exact ME reference for one chain size.

    Args:
        n_cells (int): Chain size (unit cells).

    Returns:
        tuple: list of point dicts (``ntraj``, ``mean``, ``se``, ``n``,
        ``source``) and the exact master-equation contrast (``None`` beyond
        three cells).
    """
    points = [archive_contrast(n_cells)]
    exact = None
    for path in sorted(CHECK_DIR.glob(f"solver_contrast_check_{n_cells}cells*.json")):
        block = json.loads(path.read_text())
        design = block["design"]
        if design["loss_profile"] != PROFILE or float(design["disorder_level"]) != DISORDER:
            continue
        if design.get("n_periods", 1) != 1 or design.get("kerr", 0.5) != 0.5:
            # multi-period ladder and Kerr-free stratum share the filename
            # stem but are different physics, not ladder points
            continue
        if int(design["ntraj"]) == 10:
            # production crosscheck, identical seeds and settings; the
            # archived production records already provide the N=10 point
            continue
        contrast = block["monte_carlo"]["contrast"]
        points.append({"ntraj": int(design["ntraj"]),
                       "n": int(contrast["n"]),
                       "mean": float(contrast["mean"]),
                       "se": float(contrast["se"]),
                       "source": path.name})
        me = block["exact_master_equation"]["me_contrast"]
        if me is not None:
            if exact is not None and abs(exact - float(me)) > 1e-12:
                raise ValueError(f"inconsistent ME contrast in {path}")
            exact = float(me)
    # Dedupe repeated trajectory counts (e.g. an expanded rerun sharing
    # seeds with its pilot): keep the largest sample, since overlapping
    # draws are not independent points for the fit.
    by_ntraj = {}
    for p in points:
        prev = by_ntraj.get(p["ntraj"])
        if prev is None or p["n"] > prev["n"]:
            by_ntraj[p["ntraj"]] = p
    points = sorted(by_ntraj.values(), key=lambda p: p["ntraj"])
    return points, exact


def fit_model_a(points: list[dict]) -> dict:
    """Inverse-variance weighted fit ``mean = C + c/N``.

    Args:
        points (list[dict]): Ladder points with ``ntraj``, ``mean``, ``se``.

    Returns:
        dict: ``C``, ``c``, standard error of ``C`` (None if underdetermined)
        and residual sum of squares.
    """
    x = np.asarray([1.0 / p["ntraj"] for p in points])
    y = np.asarray([p["mean"] for p in points])
    w = np.asarray([1.0 / max(p["se"], 1e-12) ** 2 for p in points])
    design = np.vstack([np.ones_like(x), x]).T
    sw = np.sqrt(w)
    (C, c), *_ = np.linalg.lstsq(design * sw[:, None], y * sw, rcond=None)
    resid = y - (C + c * x)
    dof = len(points) - 2
    se_C = None
    if dof >= 1:
        chi2 = float(np.sum(w * resid ** 2))
        cov = np.linalg.inv((design * sw[:, None]).T @ (design * sw[:, None]))
        se_C = float(np.sqrt(cov[0, 0] * chi2 / dof))
    return {"C": float(C), "c": float(c), "se_C": se_C,
            "chi2": float(np.sum(w * resid ** 2)), "dof": dof}


def fit_model_b(points: list[dict]) -> dict:
    """Unweighted fit ``mean = C + c/N + d/N^2``.

    Args:
        points (list[dict]): Ladder points with ``ntraj``, ``mean``, ``se``.

    Returns:
        dict: ``C``, ``c``, ``d`` and residual sum of squares (standard error
        of ``C`` only when the fit has at least one degree of freedom).
    """
    x = np.asarray([1.0 / p["ntraj"] for p in points])
    y = np.asarray([p["mean"] for p in points])
    design = np.vstack([np.ones_like(x), x, x ** 2]).T
    (C, c, d), *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ np.array([C, c, d])
    rss = float(np.sum(resid ** 2))
    se_C = None
    dof = len(points) - 3
    if dof >= 1:
        s2 = rss / dof
        cov = s2 * np.linalg.inv(design.T @ design)
        se_C = float(np.sqrt(cov[0, 0]))
    return {"C": float(C), "c": float(c), "d": float(d), "rss": rss,
            "se_C": se_C, "dof": dof}


def main() -> None:
    """Assemble the ladder, validate models, extrapolate, and save."""
    ladders = {}
    failures = []
    for n_cells in (2, 3, 4, 5):
        points, exact = ladder_points(n_cells)
        need = {10, 50, 200}
        have = {p["ntraj"] for p in points}
        missing = sorted(need - have)
        if missing:
            failures.append(f"{n_cells} cells: missing ntraj {missing}")
            continue
        block = {"points": points, "exact_me_contrast": exact,
                 "model_a": fit_model_a(points),
                 "model_b": fit_model_b(points)}
        if exact is not None:
            for key, model in (("model_a", block["model_a"]),
                               ("model_b", block["model_b"])):
                err = model["se_C"] or float("nan")
                z = ((model["C"] - exact) / err) if err == err else None
                model["validation"] = {"exact": exact,
                                       "C_minus_exact": model["C"] - exact,
                                       "z": z}
        ladders[f"{n_cells}cells"] = block

    # Model preference: model B must reproduce the exact ME at every
    # exactly solvable size within VALIDATION_TOL; model A is the
    # fallback if B fails anywhere.
    for key in ("model_b", "model_a"):
        bad = []
        for size, block in ladders.items():
            exact = block["exact_me_contrast"]
            if exact is None:
                continue
            model = block[key]
            if abs(model["C"] - exact) > VALIDATION_TOL:
                bad.append(f"{size}: C={model['C']:+.6f} exact={exact:+.6f}")
        if not bad:
            preferred = key
            break
        last_bad = bad
    else:
        failures.append(f"neither model validates: {last_bad}")
        preferred = None

    if failures:
        raise SystemExit("FAIL: " + "; ".join(failures))

    out = {
        "status": "SOLVER_CONTRAST_LADDER_COMPLETE",
        "preferred_model": preferred,
        "validation_tolerance": VALIDATION_TOL,
        "claim_scope": "finite-trajectory bias of the negativity estimator "
                       "and its infinite-trajectory limit; solver "
                       "convergence of the paired contrast",
        "not_supported": ["topological protection theorem",
                          "experimental validation", "metrological gain"],
        "method": {
            "N10": "archived production records (ntraj=10, seed 20260814)",
            "ladder": "run_solver_contrast_check.py at ntraj=50/200 "
                      "(same seeds, same condition)",
            "model_a": "inverse-variance weighted C + c/N",
            "model_b": "unweighted C + c/N + d/N^2",
        },
        "ladders": ladders,
    }
    OUT_PATH.write_text(json.dumps(out, indent=2) + "\n")
    for size, block in ladders.items():
        exact = block["exact_me_contrast"]
        line = f"{size}: exact ME = {exact!r}" if exact is not None else \
               f"{size}: exact ME = n/a"
        for key in ("model_a", "model_b"):
            m = block[key]
            extra = ""
            if "validation" in m:
                extra = (f"  [C-exact={m['validation']['C_minus_exact']:+.6f}"
                         f"{', z=%.2f' % m['validation']['z'] if m['validation']['z'] is not None else ''}]")
            line += f" | {key}: C={m['C']:+.6f} se={m['se_C']}{extra}"
        print(line)
    print(f"preferred model: {preferred}")
    print(f"saved {OUT_PATH}")


if __name__ == "__main__":
    main()
