#!/usr/bin/env python3
"""Adversarial-audit revision scans for the chiral Floquet production.

Extends the production protocol (run_chiral_floquet_dynamics.py) with the
scan dimensions requested by the 2026-08-17 adversarial audit without touching
the canonical production script:

  gamma      : edge-low-loss ratio scan  gamma_edge / gamma_bulk in
               {5e-4, 1e-3, 1e-2, 1e-1, 5e-1, 1.0}  (2 cells, one period)
  multiperiod: local contrast at 2T, 4T, 8T  (2 and 3 cells)
  nmax       : local-contrast panel at higher Fock cutoffs
               (2 cells N_max=3,4; 3 cells N_max=3,4)
  micromotion: observables at t = 0, T/8, ..., T (audit run R7)
  sweep      : intra/inter sweep through the gap closing (audit run R4)
  chiralloss : sublattice-localized dissipators with [S, L_k] = 0
               (audit run R5), compared against edge-low loss

Every row uses the same matched design as production: paired disorder draws,
nested trajectories, identical draw/trajectory seed scheme, edge_site even-cat
preparation, Kerr chi=0.5, 31x31 Wigner grid, ntraj=10.  The independent unit
is the disorder realization.

Checkpointing: results are written as one fragment per (loss profile, disorder
level); on resume, completed fragments are skipped.  A final merge step
assembles the fragments into the single output JSON.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from run_chiral_floquet_dynamics import (
    CANDIDATES,
    PERIOD,
    build_ops,
    driven_hamiltonian,
    initial_state,
    rate_vector,
    simulate,
)
from validate_chiral_floquet_model import evaluate_point

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SEED = 20260817
PROFILES = ("edge_low_loss", "uniform_loss")
GAMMA_RATIOS = (5e-4, 1e-3, 1e-2, 1e-1, 5e-1, 1.0)
PERIODS = (2, 4, 8)
NMAX_CELLS = {
    2: ((3, 100), (4, 100)),
    3: ((3, 100), (4, 50)),
}
NTRAJ = 10
WIGNER_GRID = 31


def fragment_dir(output: Path) -> Path:
    return output.with_suffix("") / "_fragments"


def fragment_path(fdir: Path, profile: str, disorder: float) -> Path:
    return fdir / f"{profile}_d{disorder:g}.json"


def run_config(config: dict, args: argparse.Namespace) -> None:
    n_cells = config["n_cells"]
    n_modes = 2 * n_cells
    n_max = config["n_max"]
    t_stop = config["t_stop"]
    n_time = config["n_time"]
    n_mc = config["n_mc"]
    candidates = config.get("candidates") or CANDIDATES
    profiles = tuple(config.get("profiles") or PROFILES)
    series = bool(config.get("series", False))
    ops = build_ops(n_modes, n_max)
    psi0 = initial_state(n_cells, n_max, args.initial_state_mode)
    fdir = fragment_dir(args.output)
    fdir.mkdir(parents=True, exist_ok=True)
    start = time.time()
    for profile_index, profile in enumerate(profiles):
        for disorder_index, disorder in enumerate(args.disorder_levels):
            frag = fragment_path(fdir, profile, disorder)
            if frag.is_file():
                print(f"resume: skipping {frag.name}")
                continue
            rows = []
            for sample in range(n_mc):
                draw_seed = args.seed + profile_index * 1_000_000 + disorder_index * 10_000 + sample
                rates = rate_vector(n_modes, profile, disorder, draw_seed,
                                    config["gamma_edge"], config["gamma_bulk"])
                c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
                for sector_index, (sector, candidate) in enumerate(candidates.items()):
                    seed = draw_seed + 100 + sector_index
                    row = simulate(
                        driven_hamiltonian(ops, candidate, args.kerr), psi0, c_ops, ops,
                        args.ntraj, seed, t_stop, n_time, args.wigner_grid,
                        series=series,
                    )
                    rows.append({
                        "n_cells": n_cells, "n_modes": n_modes, "n_max": n_max,
                        "sector": sector, **candidate, "kerr": args.kerr,
                        "loss_profile": profile, "disorder_level": disorder,
                        "sample": sample, "draw_seed": draw_seed,
                        "trajectory_seed": seed, "ntraj": args.ntraj,
                        "t_stop": t_stop, "period": PERIOD,
                        "wigner_grid": args.wigner_grid, "rates": rates.tolist(), **row,
                    })
            frag.write_text(json.dumps(rows, indent=2) + "\n")
            print(f"completed {frag.name} ({len(rows)} records) in {time.time() - start:.0f}s", flush=True)
    merge(args, config, start)


def merge(args: argparse.Namespace, config: dict, start: float) -> None:
    fdir = fragment_dir(args.output)
    rows = []
    for profile in tuple(config.get("profiles") or PROFILES):
        for disorder in args.disorder_levels:
            frag = fragment_path(fdir, profile, disorder)
            if not frag.is_file():
                raise FileNotFoundError(f"missing fragment {frag}")
            rows += json.loads(frag.read_text())
    output = {
        "status": "COMPLETED_REVISION_SCAN",
        "scan": args.scan,
        "claim_scope": "finite-time descriptive matched dynamics only; no protection claim",
        "not_supported": ["topological protection theorem", "experimental validation", "metrological gain"],
        "software": {"numpy": np.__version__},
        "design": {
            "n_cells": config["n_cells"], "n_max": config["n_max"],
            "initial_state_mode": args.initial_state_mode, "kerr": args.kerr,
            "n_mc": config["n_mc"], "ntraj": args.ntraj,
            "wigner_grid": args.wigner_grid,
            "gamma_edge": config["gamma_edge"], "gamma_bulk": config["gamma_bulk"],
            "t_stop": config["t_stop"], "n_time_points": config["n_time"],
            "period": PERIOD, "disorder_levels": args.disorder_levels,
            "loss_profiles": [str(p) for p in (config.get("profiles") or PROFILES)],
            "paired_disorder": True,
            "independent_unit": "disorder realization; trajectories nested",
            "seed": args.seed,
        },
        "records": rows,
        "elapsed_seconds": time.time() - start,
    }
    for key in ("budget_matched", "budget_total"):
        if key in config:
            output["design"][key] = config[key]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(f"saved {args.output} with {len(rows)} records in {output['elapsed_seconds']:.0f}s")


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="scan", required=True)

    g = sub.add_parser("gamma", help="edge/bulk loss-ratio scan")
    g.add_argument("--n-cells", type=int, default=2)
    g.add_argument("--n-mc", type=int, default=200)
    g.add_argument("--ratios", type=float, nargs="+", default=list(GAMMA_RATIOS))
    g.add_argument("--budget-matched", action="store_true",
                   help="hold the disorder-free total rate budget "
                        "sum_j gamma_j = 2*gamma_edge0 + (n_modes-2)*gamma_bulk0 "
                        "fixed across ratios (audit E14) instead of gamma_bulk")
    g.add_argument("--output-dir", type=Path, default=ROOT / "data" / "revision_scans")

    m = sub.add_parser("multiperiod", help="multi-period local contrast scan")
    m.add_argument("--n-cells", type=int, default=2)
    m.add_argument("--n-mc", type=int, default=200)
    m.add_argument("--periods", type=int, nargs="+", default=list(PERIODS))
    m.add_argument("--output-dir", type=Path, default=ROOT / "data" / "revision_scans")

    n = sub.add_parser("nmax", help="higher Fock-cutoff panel")
    n.add_argument("--n-cells", type=int, default=2)
    n.add_argument("--n-max", type=int, default=3)
    n.add_argument("--n-mc", type=int, default=100)
    n.add_argument("--output-dir", type=Path, default=ROOT / "data" / "revision_scans")

    mm = sub.add_parser("micromotion", help="sub-period observable scan (t = 0..T)")
    mm.add_argument("--n-cells", type=int, default=2)
    mm.add_argument("--n-mc", type=int, default=200)
    mm.add_argument("--output-dir", type=Path, default=ROOT / "data" / "revision_scans")

    sw = sub.add_parser("sweep", help="intra/inter sweep through the gap closing")
    sw.add_argument("--n-cells", type=int, default=2)
    sw.add_argument("--intra", type=float, nargs="+", required=True,
                    help="intra-cell hopping values of step A (units of drive frequency)")
    sw.add_argument("--total", type=float, default=1.3,
                    help="fixed intra+inter total, matching the production candidates")
    sw.add_argument("--n-mc", type=int, default=200)
    sw.add_argument("--output-dir", type=Path, default=ROOT / "data" / "revision_scans")

    cl = sub.add_parser("chiralloss", help="sublattice-localized chiral-symmetric loss")
    cl.add_argument("--n-cells", type=int, default=2)
    cl.add_argument("--n-mc", type=int, default=200)
    cl.add_argument("--output-dir", type=Path, default=ROOT / "data" / "revision_scans")

    for sp in (g, m, n, mm, sw, cl):
        sp.add_argument("--disorder-levels", type=float, nargs="+", default=(0.0, 0.5))
        sp.add_argument("--ntraj", type=int, default=NTRAJ)
        sp.add_argument("--wigner-grid", type=int, default=WIGNER_GRID)
        sp.add_argument("--initial-state-mode", choices=("edge_site", "topological_zero_mode"), default="edge_site")
        sp.add_argument("--kerr", type=float, default=0.5)
        sp.add_argument("--gamma-edge", type=float, default=1e-4)
        sp.add_argument("--gamma-bulk", type=float, default=0.2)
        sp.add_argument("--seed", type=int, default=DEFAULT_SEED)
        sp.add_argument("--t-stop", type=float, default=PERIOD)
        sp.add_argument("--profiles", nargs="+", default=None,
                        help="loss profiles (default: edge_low_loss uniform_loss)")
        sp.add_argument("--output", type=Path, default=None)

    args = p.parse_args()

    configs = []
    if args.scan == "gamma":
        n_modes = 2 * args.n_cells
        if args.budget_matched:
            # Disorder-free base-rate sum of the production edge_low_loss
            # profile: two edge sites at gamma_edge0, the rest at gamma_bulk0.
            budget = 2.0 * args.gamma_edge + (n_modes - 2) * args.gamma_bulk
        for ratio in args.ratios:
            if args.budget_matched:
                # gamma_edge = ratio * gamma_bulk with
                # (n_modes-2)*gamma_bulk + 2*gamma_edge = budget exactly.
                gamma_bulk = budget / ((n_modes - 2) + 2.0 * ratio)
                gamma_edge = ratio * gamma_bulk
                kind = f"matched_{ratio:g}"
            else:
                gamma_bulk = args.gamma_bulk
                gamma_edge = args.gamma_bulk * ratio
                kind = f"gamma_ratio_{ratio:g}"
            config = {
                "kind": kind,
                "n_cells": args.n_cells, "n_max": 2,
                "gamma_edge": gamma_edge, "gamma_bulk": gamma_bulk,
                "t_stop": PERIOD, "n_time": 9, "n_mc": args.n_mc,
            }
            if args.budget_matched:
                config["budget_matched"] = True
                config["budget_total"] = budget
            configs.append(config)
    elif args.scan == "multiperiod":
        for np_ in args.periods:
            configs.append({
                "kind": f"multiperiod_{np_}T",
                "n_cells": args.n_cells, "n_max": 2,
                "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
                "t_stop": np_ * PERIOD, "n_time": 1 + 8 * np_, "n_mc": args.n_mc,
            })
    elif args.scan == "nmax":
        configs.append({
            "kind": f"nmax{args.n_max}",
            "n_cells": args.n_cells, "n_max": args.n_max,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": PERIOD, "n_time": 9, "n_mc": args.n_mc,
        })
    elif args.scan == "micromotion":
        configs.append({
            "kind": "micromotion",
            "n_cells": args.n_cells, "n_max": 2,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": PERIOD, "n_time": 9, "n_mc": args.n_mc, "series": True,
        })
    elif args.scan == "sweep":
        for intra in args.intra:
            inter = args.total - intra
            if inter <= 0:
                p.error("each intra value must be smaller than --total")
            audit = evaluate_point(intra, inter)
            last = audit["resolution_results"][-1]
            candidate = {"sweep_point": {
                "intra_A": intra, "inter_B": inter,
                "drive_ratio": intra / inter,
                "nu_0": last["nu_0"], "nu_pi": last["nu_pi"],
                "winding_frame_A": last["frame_A"]["winding_integer"],
                "winding_frame_B": last["frame_B"]["winding_integer"],
                "resolution_stable": audit["resolution_stable"],
                "bulk_point_candidate_pass": audit["bulk_point_candidate_pass"],
                "minimum_zero_gap": audit["floquet_gaps"]["minimum_zero_gap"],
                "minimum_pi_gap": audit["floquet_gaps"]["minimum_pi_gap"],
            }}
            configs.append({
                "kind": f"sweep_intra{intra:g}",
                "n_cells": args.n_cells, "n_max": 2,
                "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
                "t_stop": PERIOD, "n_time": 9, "n_mc": args.n_mc,
                "candidates": candidate, "profiles": ("edge_low_loss",),
            })
    elif args.scan == "chiralloss":
        configs.append({
            "kind": "chiralloss",
            "n_cells": args.n_cells, "n_max": 2,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": PERIOD, "n_time": 9, "n_mc": args.n_mc,
            "profiles": tuple(args.profiles or
                              ("edge_low_loss", "chiral_even_loss",
                               "chiral_odd_loss", "chiral_even_loss_matched")),
        })

    for cfg in configs:
        if args.profiles is not None and "profiles" not in cfg:
            cfg["profiles"] = tuple(args.profiles)
        if args.output is None:
            args.output = args.output_dir / f"{args.scan}_{cfg['kind']}_{cfg['n_cells']}cells.json"
        args.output.parent.mkdir(parents=True, exist_ok=True)
        run_config(cfg, args)
        args.output = None


if __name__ == "__main__":
    main()
