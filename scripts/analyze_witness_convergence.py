#!/usr/bin/env python3
"""Witness alpha-search convergence (audit follow-up: alpha-search study).

The production witness analysis (``scripts/analyze_r8_witness.py``)
maximises each Fock fidelity ``F_n(alpha)`` over an 81-point coarse grid
on ``[-2, 2]^2`` plus a local refinement (step 0.0025) around the single
coarse argmax.  This script re-maximises on the identical dumped 3x3
blocks with a denser search -- a 161-point global grid (step 0.025) plus
a 0.0005-step refinement around the top-3 coarse maxima -- and reports,
per stratum, the shift of the witness maxima, certification-count flips,
and the paired ``Delta F_1`` shift.

Gates (script aborts on failure):
  1. every record of the three rerun files is processed;
  2. the dense maximum is never below the coarse-81 maximum (the coarse
     grid is a subset of the dense grid, so this is a mathematical
     identity and gates the implementation).

Output:
    data/witness_alpha_convergence.json
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import time
from pathlib import Path

import numpy as np

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_r8_witness import (  # noqa: E402
    N_WITNESS,
    OMEGA,
    block_matrix,
    maximize_witness,
    wilson_interval,
    witness_fidelity,
)

ROOT = Path(__file__).resolve().parents[1]
DENSE = np.linspace(-2.0, 2.0, 161)
REFINE_STEP = 0.0005
REFINE_HALF_WIDTH = 0.05
TOP_K = 3

_RECORDS: list = []


def maximize_dense(rho: np.ndarray, n: int) -> tuple[float, complex]:
    """Dense alpha maximisation of ``F_n`` with top-3 local refinement.

    Args:
        rho (ndarray): 3x3 reduced density matrix (dimensionless).
        n (int): Witness Fock index (dimensionless).

    Returns:
        tuple: ``(max F_n, argmax alpha)`` (dimensionless).
    """
    re, im = np.meshgrid(DENSE, DENSE, indexing="ij")
    a0 = (re + 1j * im).ravel()
    f0 = witness_fidelity(rho, n, a0)
    order = np.argsort(f0)[::-1][:TOP_K]
    best = float(f0[order[0]])
    best_a = complex(a0[order[0]])
    step = DENSE[1] - DENSE[0]
    for i in order:
        c = complex(a0[i])
        gr = np.arange(c.real - REFINE_HALF_WIDTH,
                       c.real + REFINE_HALF_WIDTH + REFINE_STEP, REFINE_STEP)
        gi = np.arange(c.imag - REFINE_HALF_WIDTH,
                       c.imag + REFINE_HALF_WIDTH + REFINE_STEP, REFINE_STEP)
        rr, ii = np.meshgrid(gr, gi, indexing="ij")
        a1 = (rr + 1j * ii).ravel()
        f1 = witness_fidelity(rho, n, a1)
        j = int(np.argmax(f1))
        if float(f1[j]) > best:
            best, best_a = float(f1[j]), complex(a1[j])
    _ = step
    return best, best_a


def maximize_coarse81(rho: np.ndarray, n: int) -> float:
    """Coarse-81 global maximum (no refinement): the subset bound.

    Args:
        rho (ndarray): 3x3 reduced density matrix (dimensionless).
        n (int): Witness Fock index (dimensionless).

    Returns:
        float: max of ``F_n`` over the original 81-point grid.
    """
    g = np.linspace(-2.0, 2.0, 81)
    re, im = np.meshgrid(g, g, indexing="ij")
    return float(witness_fidelity(rho, n, (re + 1j * im).ravel()).max())


def _init_worker(records: list) -> None:
    """Share the record list with forked workers.

    Args:
        records (list): Loaded rerun records with parsed blocks.
    """
    global _RECORDS
    _RECORDS = records


def _task(idx: int) -> tuple:
    """Dense + production + coarse maxima for one record, both sites, all n.

    Args:
        idx (int): Index into the shared record list.

    Returns:
        tuple: ``(idx, {site: {n: (dense, prod, coarse81)}})``.
    """
    r = _RECORDS[idx]
    out = {}
    for site_field in ("rho_edge", "rho_center"):
        rho = block_matrix(r[site_field])
        per_n = {}
        for n in N_WITNESS:
            dense, _ = maximize_dense(rho, n)
            prod, _ = maximize_witness(rho, n)
            coarse = maximize_coarse81(rho, n)
            per_n[n] = (dense, prod, coarse)
        out[site_field] = per_n
    return idx, out


def main() -> None:
    """Run the dense search over the full rerun ensemble and aggregate."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--output", type=Path,
                   default=ROOT / "data" / "witness_alpha_convergence.json")
    args = p.parse_args()

    start = time.time()
    records = []
    for cells in (2, 3, 4):
        rerun = json.loads((ROOT / "data" / "r8_witness"
                            / f"r8_witness_{cells}cells.json").read_text())
        for r in rerun["records"]:
            r["_cells"] = cells
            records.append(r)
    print(f"loaded {len(records)} records", flush=True)

    ctx = multiprocessing.get_context("fork")
    results = {}
    with ctx.Pool(args.workers, initializer=_init_worker,
                  initargs=(records,)) as pool:
        for i, (idx, res) in enumerate(pool.imap_unordered(_task, range(len(records)))):
            results[idx] = res
            if (i + 1) % 2000 == 0:
                print(f"  {i + 1}/{len(records)}", flush=True)

    gates = {"all_processed": len(results) == len(records),
             "dense_ge_coarse": True}
    worse = 0
    prod_flips = 0
    max_prod_delta = 0.0
    for idx, res in results.items():
        for site_field in ("rho_edge", "rho_center"):
            for n in N_WITNESS:
                dense, prod, coarse = res[site_field][n]
                if not (np.isfinite(dense) and np.isfinite(prod)
                        and np.isfinite(coarse)):
                    gates["dense_ge_coarse"] = False
                if dense < coarse - 1e-12:
                    worse += 1
                if abs(dense - prod) > 1e-4:
                    prod_flips += 1
                max_prod_delta = max(max_prod_delta, abs(dense - prod))
    if worse:
        gates["dense_ge_coarse"] = False

    summary = {}
    max_shift = 0.0
    for cells in (2, 3, 4):
        for disorder in sorted({r["disorder_level"] for r in records
                                if r["_cells"] == cells}):
            for profile in sorted({r["loss_profile"] for r in records
                                   if r["_cells"] == cells}):
                sel = [i for i, r in enumerate(records)
                       if (r["_cells"], r["disorder_level"],
                           r["loss_profile"]) == (cells, disorder, profile)]
                entry = {"n": len(sel) // 2}
                for sector in ("chiral_topological", "chiral_trivial"):
                    grp = [i for i in sel
                           if records[i]["sector"] == sector]
                    vals_o = np.array([results[i]["rho_edge"][1][1]
                                       for i in grp])
                    vals_d = np.array([results[i]["rho_edge"][1][0]
                                       for i in grp])
                    lo, hi = OMEGA[1]
                    co = int(np.sum(vals_o > hi))
                    cd = int(np.sum(vals_d > hi))
                    entry[sector] = {
                        "certified_orig_proxy": co,
                        "certified_dense": cd,
                        "certified_flips": int(cd - co),
                        "f1_mean_orig_proxy": float(vals_o.mean()),
                        "f1_mean_dense": float(vals_d.mean()),
                        "f1_shift": float(vals_d.mean() - vals_o.mean()),
                        "f1_max_shift": float(np.abs(vals_d - vals_o).max()),
                    }
                    max_shift = max(max_shift,
                                    float(np.abs(vals_d - vals_o).max()))
                topo = [i for i in sel
                        if records[i]["sector"] == "chiral_topological"]
                triv = [i for i in sel
                        if records[i]["sector"] == "chiral_trivial"]
                do = np.array([results[t]["rho_edge"][1][1]
                               - results[v]["rho_edge"][1][1]
                               for t, v in zip(topo, triv)])
                dd = np.array([results[t]["rho_edge"][1][0]
                               - results[v]["rho_edge"][1][0]
                               for t, v in zip(topo, triv)])
                entry["paired_dF1"] = {
                    "mean_orig_proxy": float(do.mean()),
                    "mean_dense": float(dd.mean()),
                    "shift": float(dd.mean() - do.mean()),
                }
                summary[f"{cells}cells_d{disorder:g}_{profile}"] = entry

    # centre-site control (n=1..6): does any dense maximum cross threshold?
    center_cross = {}
    for n in N_WITNESS:
        lo, hi = OMEGA[n]
        mx = max(results[i]["rho_center"][n][0] for i in results)
        center_cross[f"n{n}"] = {"dense_global_max": float(mx),
                                 "omega_hi": hi,
                                 "crosses": bool(mx > hi)}

    out = {
        "status": ("WITNESS_ALPHA_CONVERGENCE_COMPLETE"
                   if all(gates.values())
                   else "WITNESS_ALPHA_CONVERGENCE_FAILED"),
        "claim_scope": "alpha-grid convergence of the grid-free witness "
                       "maxima on the identical dumped blocks; search "
                       "validation only",
        "not_supported": ["topological protection theorem",
                          "experimental validation", "metrological gain"],
        "design": {
            "coarse_grid": "81 points on [-2,2], step 0.05",
            "dense_grid": "161 points on [-2,2], step 0.025",
            "refinement": f"top-{TOP_K} maxima, +- {REFINE_HALF_WIDTH} "
                          f"at step {REFINE_STEP}",
            "n_records": len(records),
        },
        "gates": gates,
        "worse_than_coarse_count": worse,
        "max_abs_f1_shift": max_shift,
        "prod_vs_dense_gt_1e-4_count": prod_flips,
        "max_abs_f1_delta_vs_prod": max_prod_delta,
        "strata": summary,
        "center_control": center_cross,
        "software": {"numpy": np.__version__},
        "elapsed_seconds": time.time() - start,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(f"gates: {gates}; worse={worse}; max_shift={max_shift:.3e}")
    print(f"saved {args.output}")
    if not all(gates.values()):
        raise SystemExit("FAIL: " + ", ".join(k for k, v in gates.items()
                                              if not v))


if __name__ == "__main__":
    main()
