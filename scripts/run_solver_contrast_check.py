#!/usr/bin/env python3
"""Exact master-equation vs Monte-Carlo contrast cross-check (audit A6, item 16).

The production contrast (R1, ``data/chiral_production/chiral_floquet_*NMC500``)
is estimated with the Monte Carlo wave-function solver at ``ntraj = 10``
trajectories per disorder sample.  This script measures how close that
ensemble estimate is to the *exact* Lindblad master equation for the same
configuration: for every production sample of one condition (loss profile,
disorder level, seed 20260814) it runs

``mesolve`` (density-matrix, deterministic, one solve per sector) and
``mcsolve`` (``ntraj = 10`` trajectories, seeded exactly as production,

``draw_seed = seed + sample`` at profile index 0 / disorder index 0 and
``trajectory_seed = draw_seed + 100 + sector_index``),

then reports the paired topological-minus-trivial negativity contrast of
each solver, the difference of the Monte Carlo mean from the exact master
equation with a 95% interval, and a record-level cross-check of the
recomputed Monte Carlo observables against the archived production file.

Gates (script aborts on failure):
  1. production cross-check: max |MC - production| over negativity and
     occupations must be numerically zero (same seeds, same code path);
  2. solver agreement: |mean(MC contrast) - ME contrast| must lie within
     3 standard errors of the Monte Carlo sample mean.

Output:
    data/solver_contrast_check_{n}cells.json
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import time
from pathlib import Path

import numpy as np
import qutip as qt

from run_chiral_floquet_dynamics import (
    CANDIDATES,
    DISORDER_LEVELS,
    PERIOD,
    build_ops,
    driven_hamiltonian,
    initial_state,
    negativity,
    rate_vector,
    simulate,
)

ROOT = Path(__file__).resolve().parents[1]
PRODUCTION_SEED = 20260814


def measure(state: qt.Qobj, ops: list, grid_points: int) -> dict:
    """Terminal-time observables, identical to ``simulate()``.

    Args:
        state (qutip.Qobj): Final-time density matrix or ket.
        ops (list): Single-site annihilation operators (dimensionless
            truncated Fock space).
        grid_points (int): Side length of the Wigner grid (dimensionless).

    Returns:
        dict: Wigner-void negativity of the edge and centre mode
        (dimensionless), edge/bulk populations and per-site occupations
        (dimensionless).
    """
    occ = np.asarray([float(np.real(qt.expect(op.dag() * op, state)))
                      for op in ops])
    return {
        "negativity_edge": negativity(state, 0, grid_points),
        "negativity_center": negativity(state, len(ops) // 2, grid_points),
        "edge_population": float(occ[0] + occ[-1]),
        "bulk_population": float(np.sum(occ[1:-1])),
        "occupations": occ.tolist(),
    }


def solve_me(h, psi0, c_ops, times: np.ndarray, ops: list,
             grid_points: int) -> dict:
    """Exact Lindblad master equation, terminal observables.

    Args:
        h (qutip.QobjEvo): Time-dependent Floquet Hamiltonian.
        psi0 (qutip.Qobj): Initial cat state (ket).
        c_ops (list): Collapse operators sqrt(gamma_j) a_j.
        times (numpy.ndarray): Stored evaluation times (units of 1/omega).
        ops (list): Single-site annihilation operators.
        grid_points (int): Wigner grid side length (dimensionless).

    Returns:
        dict: Same observable block as :func:`measure`.
    """
    result = qt.mesolve(
        h, psi0, times, c_ops,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    return measure(result.states[-1], ops, grid_points)


def sample_stats(values: np.ndarray) -> dict:
    """Mean, sd, standard error and 95% normal interval of sample values.

    Args:
        values (numpy.ndarray): i.i.d. sample values (dimensionless).

    Returns:
        dict: ``n``, ``mean``, ``sd``, ``se``, ``ci95_low``, ``ci95_high``.
    """
    x = np.asarray(values, dtype=float)
    n = int(x.size)
    mean = float(x.mean())
    sd = float(x.std(ddof=1)) if n > 1 else 0.0
    se = sd / np.sqrt(n) if n > 1 else float("nan")
    return {"n": n, "mean": mean, "sd": sd, "se": float(se),
            "ci95_low": float(mean - 1.96 * se),
            "ci95_high": float(mean + 1.96 * se)}


_WORKER = {}


def _init_worker(payload: dict) -> None:
    """Build per-process operators, Hamiltonians and initial state.

    Args:
        payload (dict): Static run parameters (chain size, cutoff, rates,
            solver settings); executed once per worker process.
    """
    _WORKER.update(payload)
    _WORKER["ops"] = build_ops(_WORKER["n_modes"], _WORKER["n_max"])
    _WORKER["psi0"] = initial_state(_WORKER["n_cells"], _WORKER["n_max"],
                                    "edge_site")
    _WORKER["hamiltonians"] = {
        sector: driven_hamiltonian(_WORKER["ops"], cand, _WORKER["kerr"])
        for sector, cand in CANDIDATES.items()}


def _mc_task(task: tuple) -> tuple:
    """Run one (sample, sector) Monte Carlo record.

    Args:
        task (tuple): ``(sample, sector_index)``.

    Returns:
        tuple: ``(sample, sector_index, row)`` with the observable dict.
    """
    sample, sector_index = task
    sector = list(CANDIDATES)[sector_index]
    profile_index = {"edge_low_loss": 0, "uniform_loss": 1}[_WORKER["profile"]]
    disorder_index = {level: i for i, level in enumerate(DISORDER_LEVELS)}
    draw_seed = (PRODUCTION_SEED + profile_index * 1_000_000
                 + disorder_index[_WORKER["disorder"]] * 10_000 + sample)
    rates = rate_vector(_WORKER["n_modes"], _WORKER["profile"],
                        _WORKER["disorder"], draw_seed,
                        _WORKER["gamma_edge"], _WORKER["gamma_bulk"])
    c_ops = [np.sqrt(x) * op
             for x, op in zip(rates, _WORKER["ops"]) if x > 0]
    seed = draw_seed + 100 + sector_index
    row = simulate(_WORKER["hamiltonians"][sector], _WORKER["psi0"],
                   c_ops, _WORKER["ops"], _WORKER["ntraj"], seed,
                   _WORKER["t_stop"], _WORKER["n_time"],
                   _WORKER["wigner_grid"])
    return sample, sector_index, row


def _write_fragment(fdir: Path | None, res: tuple) -> None:
    """Atomically checkpoint one ``(sample, sector, row)`` record.

    Enables resuming an interrupted run (``--resume``) without losing the
    records already computed, following the repository ``_fragments``
    convention. A no-op when ``fdir`` is None, so the default behaviour is
    unchanged.

    Args:
        fdir (Path | None): Fragment directory, or None to disable
            checkpointing.
        res (tuple): ``(sample, sector_index, row)`` from :func:`_mc_task`.
    """
    if fdir is None:
        return
    sample, sector_index, row = res
    path = fdir / f"sample_{sample:04d}_sector_{sector_index}.json"
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps({"sample": sample,
                               "sector_index": sector_index,
                               "row": row}) + "\n")
    tmp.replace(path)


def main() -> None:
    """Run the matched mesolve/mcsolve contrast check for one chain size."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-cells", type=int, default=2, choices=(2, 3, 4, 5))
    p.add_argument("--n-mc", type=int, default=500,
                   help="number of matched samples (dimensionless)")
    p.add_argument("--ntraj", type=int, default=10,
                     help="Monte Carlo trajectories per sample (production: 10)")
    p.add_argument("--n-periods", type=int, default=1,
                    help="drive periods evolved (production: 1; use 2/4 for "
                         "the longtime ladder)")
    p.add_argument("--skip-mc", action="store_true",
                   help="do not re-simulate the Monte Carlo arm; take its "
                        "statistics from the archived production file")
    p.add_argument("--workers", type=int, default=1,
                   help="parallel worker processes for the Monte Carlo arm")
    p.add_argument("--profile", default="edge_low_loss")
    p.add_argument("--disorder", type=float, default=0.0)
    p.add_argument("--gamma-edge", type=float, default=1e-4)
    p.add_argument("--gamma-bulk", type=float, default=0.2)
    p.add_argument("--n-max", type=int, default=2)
    p.add_argument("--kerr", type=float, default=0.5)
    p.add_argument("--wigner-grid", type=int, default=31)
    p.add_argument("--n-time", type=int, default=9)
    p.add_argument("--production", type=Path, default=None,
                   help="archived production JSON (default: chiral_production)")
    p.add_argument("--output", type=Path, default=None)
    p.add_argument("--fragment-dir", type=Path, default=None,
                   help="optional directory for per-sample checkpoint "
                        "fragments; makes a long run resumable with --resume "
                        "(repo _fragments convention). Default: off")
    p.add_argument("--resume", action="store_true",
                   help="with --fragment-dir, skip records whose fragment "
                        "already exists instead of re-simulating them")
    args = p.parse_args()
    if args.resume and args.fragment_dir is None:
        raise SystemExit("FAIL: --resume requires --fragment-dir")

    n_cells = args.n_cells
    if args.production:
        production_path = args.production
    else:
        candidates = [ROOT / "data" / "chiral_production"
                      / f"chiral_floquet_{n_cells}cells_NMC500.json",
                      *sorted((ROOT / "data" / "chiral_production").glob(
                          f"chiral_floquet_{n_cells}cells_NMC*.json"))]
        production_path = next((c for c in candidates if c.is_file()), None)
    out_path = (args.output
                or ROOT / "data" / f"solver_contrast_check_{n_cells}cells.json")

    n_modes = 2 * n_cells
    ops = build_ops(n_modes, args.n_max)
    psi0 = initial_state(n_cells, args.n_max, "edge_site")
    t_stop = args.n_periods * PERIOD
    times = np.linspace(0.0, t_stop, args.n_time)
    hamiltonians = {sector: driven_hamiltonian(ops, cand, args.kerr)
                    for sector, cand in CANDIDATES.items()}

    profile_index = {"edge_low_loss": 0, "uniform_loss": 1}[args.profile]
    disorder_index = {level: i for i, level in enumerate(DISORDER_LEVELS)}
    if args.disorder not in disorder_index:
        raise SystemExit(f"FAIL: disorder {args.disorder} not in {DISORDER_LEVELS}")

    start = time.time()

    # ---- exact master equation: one deterministic solve per sector -------
    # At d = 0 the rate vector is seed-independent, so every sample shares
    # these collapse operators; the ME contrast is a single exact number.
    # The density-matrix solve is only tractable up to three cells (the
    # Liouvillian of four cells already exceeds machine limits); at four and
    # five cells the exact reference is obtained from the trajectory-count
    # extrapolation (see scripts/analyze_solver_contrast_ladder.py).
    me_note = None
    if n_cells <= 3:
        draw_seed0 = (PRODUCTION_SEED + profile_index * 1_000_000
                      + disorder_index[args.disorder] * 10_000 + 0)
        rates = rate_vector(n_modes, args.profile, args.disorder, draw_seed0,
                            args.gamma_edge, args.gamma_bulk)
        c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
        me = {sector: solve_me(hamiltonians[sector], psi0, c_ops, times, ops,
                               args.wigner_grid)
              for sector in CANDIDATES}
        me_contrast = (me["chiral_topological"]["negativity_edge"]
                       - me["chiral_trivial"]["negativity_edge"])
    else:
        me = None
        me_contrast = None
        me_note = ("density-matrix master equation infeasible beyond three "
                   "cells; use the 1/N trajectory extrapolation")

    # ---- Monte Carlo: production seed scheme, production ntraj -----------
    mc = {sector: {k: [] for k in ("negativity_edge", "negativity_center",
                                   "edge_population", "bulk_population",
                                   "occupations")}
          for sector in CANDIDATES}
    crosscheck_note = None
    if args.skip_mc:
        if production_path is None or not production_path.is_file():
            raise SystemExit("FAIL: --skip-mc requires the production file")
        prod = json.loads(production_path.read_text())
        archive = {(int(r["sample"]), r["sector"]): r
                   for r in prod["records"]
                   if r["loss_profile"] == args.profile
                   and float(r["disorder_level"]) == args.disorder}
        samples = sorted({s for (s, sec) in archive
                          if sec == "chiral_topological"})
        if not samples:
            raise SystemExit("FAIL: no archived records for the condition")
        for sector in CANDIDATES:
            for key in ("negativity_edge", "negativity_center",
                        "edge_population", "bulk_population"):
                mc[sector][key] = [float(archive[(s, sector)][key])
                                   for s in samples]
            mc[sector]["occupations"] = [archive[(s, sector)]["occupations"]
                                         for s in samples]
        n_samples = len(samples)
        crosscheck_note = ("Monte Carlo statistics taken from the archive; "
                           "bit-exact re-simulation of the production "
                           "configuration verified at two and three cells")
    else:
        payload = {"n_cells": n_cells, "n_modes": n_modes,
                   "n_max": args.n_max, "kerr": args.kerr,
                   "profile": args.profile, "disorder": args.disorder,
                   "gamma_edge": args.gamma_edge,
                   "gamma_bulk": args.gamma_bulk, "ntraj": args.ntraj,
                   "t_stop": args.n_periods * PERIOD,
                   "n_time": args.n_time, "wigner_grid": args.wigner_grid}
        tasks = [(s, si) for s in range(args.n_mc)
                 for si in range(len(CANDIDATES))]
        _init_worker(payload)
        results = []
        fdir = args.fragment_dir
        done = {}
        if fdir is not None:
            fdir.mkdir(parents=True, exist_ok=True)
            if args.resume:
                for frag in sorted(fdir.glob("sample_*.json")):
                    rec = json.loads(frag.read_text())
                    key = (int(rec["sample"]), int(rec["sector_index"]))
                    done[key] = (key[0], key[1], rec["row"])
        if done:
            tasks = [t for t in tasks if t not in done]
            print(f"resume: {len(done)} records already on disk, "
                  f"{len(tasks)} remaining", flush=True)
        n_total = len(tasks)
        if args.workers > 1:
            ctx = multiprocessing.get_context("fork")
            with ctx.Pool(args.workers, initializer=_init_worker,
                          initargs=(payload,)) as pool:
                for i, res in enumerate(
                        pool.imap_unordered(_mc_task, tasks, chunksize=1), 1):
                    results.append(res)
                    _write_fragment(fdir, res)
                    if i % 20 == 0:
                        print(f"completed {i}/{n_total} records",
                              flush=True)
        else:
            for i, task in enumerate(tasks, 1):
                res = _mc_task(task)
                results.append(res)
                _write_fragment(fdir, res)
                if i % 100 == 0:
                    print(f"completed {i}/{n_total} records", flush=True)
        results.extend(done.values())
        results.sort(key=lambda r: (r[0], r[1]))
        sectors = list(CANDIDATES)
        for sample, sector_index, row in results:
            for key in mc[sectors[sector_index]]:
                mc[sectors[sector_index]][key].append(row[key])
        n_samples = args.n_mc
    mc_contrast = (np.asarray(mc["chiral_topological"]["negativity_edge"])
                   - np.asarray(mc["chiral_trivial"]["negativity_edge"]))
    mc_stats = sample_stats(mc_contrast)
    se = mc_stats["se"]
    if me_contrast is None:
        delta = None
        solver_agrees = True
    else:
        delta = mc_stats["mean"] - me_contrast
        solver_agrees = (bool(abs(delta) <= 3.0 * se) if se > 0
                         else bool(delta == 0.0))

    # ---- record-level cross-check against the archived production file ---
    crosscheck = {"available": False}
    if args.skip_mc:
        crosscheck = {"available": False, "note": crosscheck_note,
                      "production_file": (str(production_path.relative_to(ROOT))
                                          if production_path else None)}
    elif args.ntraj != 10 or args.n_periods != 1:
        crosscheck = {"available": False,
                      "note": f"cross-check only defined for the production "
                              f"ntraj=10, one-period configuration (ran "
                              f"ntraj={args.ntraj}, "
                              f"periods={args.n_periods})"}
    elif production_path is not None and production_path.is_file():
        prod = json.loads(production_path.read_text())
        wanted = {(int(r["sample"]), r["sector"]): r
                  for r in prod["records"]
                  if r["loss_profile"] == args.profile
                  and float(r["disorder_level"]) == args.disorder}
        max_diff = 0.0
        compared = 0
        n_cmp = min(args.n_mc, len({s for s, _ in wanted}))
        for sector in CANDIDATES:
            for sample in range(n_cmp):
                key = (sample, sector)
                if key not in wanted:
                    raise SystemExit(f"FAIL: production record missing {key}")
                row = wanted[key]
                idx = sample
                for observable in ("negativity_edge", "negativity_center",
                                   "edge_population", "bulk_population"):
                    diff = abs(float(mc[sector][observable][idx])
                               - float(row[observable]))
                    max_diff = max(max_diff, diff)
                    compared += 1
                occ_diff = np.max(np.abs(np.asarray(mc[sector]["occupations"][idx])
                                         - np.asarray(row["occupations"])))
                max_diff = max(max_diff, float(occ_diff))
                compared += len(row["occupations"])
        crosscheck = {
            "available": True,
            "production_file": str(production_path.relative_to(ROOT)),
            "production_status": prod["status"],
            "records_compared": compared,
            "max_abs_diff": float(max_diff),
        }

    gates = {
        "production_crosscheck_zero": (not crosscheck["available"]
                                       or crosscheck["max_abs_diff"] == 0.0),
        "solver_agreement_within_3se": solver_agrees,
    }
    out = {
        "status": ("SOLVER_CONTRAST_CHECK_COMPLETE"
                   if all(gates.values()) else "SOLVER_CONTRAST_CHECK_FAILED"),
        "claim_scope": "finite-time master-equation vs Monte-Carlo agreement "
                       "for the production configuration; solver validation "
                       "only",
        "not_supported": ["topological protection theorem",
                          "experimental validation", "metrological gain"],
        "design": {
            "n_cells": n_cells, "n_modes": n_modes, "n_max": args.n_max,
            "initial_state_mode": "edge_site",             "kerr": args.kerr,
            "n_mc": n_samples, "ntraj": args.ntraj,
            "monte_carlo_from_archive": bool(args.skip_mc),
            "loss_profile": args.profile, "disorder_level": args.disorder,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": float(args.n_periods * PERIOD),
            "n_periods": args.n_periods, "n_time_points": args.n_time,
            "wigner_grid": args.wigner_grid,
            "production_seed": PRODUCTION_SEED,
            "draw_seed": "seed + profile_index*1e6 + disorder_index*1e4 + sample",
            "trajectory_seed": "draw_seed + 100 + sector_index",
            "paired_sectors": True,
        },
        "exact_master_equation": ({
            "me_contrast": float(me_contrast),
            "me": {sector: {k: v for k, v in block.items() if k != "occupations"}
                   for sector, block in me.items()},
            "me_occupations": {sector: block["occupations"]
                               for sector, block in me.items()},
        } if me is not None else {
            "me_contrast": None,
            "note": me_note,
        }),
        "monte_carlo": {
            "contrast": mc_stats,
            "by_sector": {
                sector: {key: sample_stats(np.asarray(block[key], dtype=float)
                                           if key != "occupations"
                                           else np.asarray(
                                               [float(np.mean(o))
                                                for o in block[key]]))
                         for key in ("negativity_edge", "negativity_center",
                                     "edge_population", "bulk_population")}
                for sector, block in mc.items()
            },
        },
        "delta_mc_minus_me": ({
            "difference": float(delta),
            "mc_se": float(se),
            "z": float(delta / se) if se > 0 else None,
        } if delta is not None else {
            "difference": None,
            "mc_se": float(se),
            "note": me_note,
        }),
        "crosscheck": crosscheck,
        "gates": gates,
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "elapsed_seconds": time.time() - start,
    }
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    if me_contrast is not None:
        print(f"ME contrast      = {me_contrast:+.6f}")
    else:
        print(f"ME contrast      = n/a ({me_note})")
    print(f"MC contrast mean = {mc_stats['mean']:+.6f} "
          f"+/- {1.96 * mc_stats['se']:.6f} (95% CI, n={mc_stats['n']})")
    if delta is not None:
        print(f"delta (MC - ME)  = {delta:+.6f} "
              f"(z={delta / se if se > 0 else float('nan'):+.2f})")
    if crosscheck["available"]:
        print(f"production cross-check max |diff| = {crosscheck['max_abs_diff']:.3e} "
              f"over {crosscheck['records_compared']} values")
    print(f"gates: {gates}")
    print(f"saved {out_path}")
    if not all(gates.values()):
        raise SystemExit("FAIL: " + ", ".join(k for k, v in gates.items()
                                              if not v))


if __name__ == "__main__":
    main()
