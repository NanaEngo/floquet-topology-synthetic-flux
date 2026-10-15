#!/usr/bin/env python3
"""Loss-ratio x size matrix for the Track A sign pattern (A-S2).

Reuses run_collective_wigner.run_dynamics (sector-matched chiral
candidates, paired disorder draws) over loss ratios
gamma_edge/gamma_bulk in {5e-4, 1e-2, 0.1, 0.5} at fixed bulk and cell
counts {2} (+3 with --extended), recording the paired
topological-minus-trivial site-0 contrast sign per cell. Locates the
sign-flip boundary between the production ratio and the loss-dominated
regime. HPC-scale at full sizes/samples; --smoke runs one cell.
Output: data/ratio_size_matrix.json.
"""
from __future__ import annotations

import argparse
import json
import sys
import types
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_collective_wigner as rcw

OUT = ROOT / "data" / "ratio_size_matrix.json"
RATIOS = (5e-4, 1e-2, 0.1, 0.5)
GAMMA_BULK = 0.2


def paired_contrast(records):
    """Mean topological-minus-trivial site-0 contrast over matched draws.

    Args:
        records: Collective-wigner record list with sector, sample,
            disorder_level and negativity_local_site0 fields.

    Returns:
        Tuple (mean, se, n) in Wigner-negativity units (dimensionless).
    """
    groups = {}
    for r in records:
        groups.setdefault((r["sample"], r["disorder_level"]), {})[r["sector"]] = r
    vals = [g["chiral_topological"]["negativity_local_site0"]
            - g["chiral_trivial"]["negativity_local_site0"]
            for g in groups.values() if len(g) == 2]
    vals = np.asarray(vals)
    n = len(vals)
    return float(vals.mean()), float(vals.std(ddof=1) / n ** 0.5) if n > 1 else 0.0, n


def main() -> int:
    """Run the ratio x size matrix with per-cell checkpointing."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--extended", action="store_true",
                    help="include three cells (HPC-scale)")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    cells = (2,) if (args.smoke or not args.extended) else (2, 3)
    n_mc, ntraj = (2, 4) if args.smoke else (10, 10)
    try:
        matrix = json.loads(args.out.read_text())["cells"]
        done = {(c["n_cells"], c["ratio"]) for c in matrix}
    except (FileNotFoundError, KeyError):
        matrix, done = [], set()
    for n_cells in cells:
        for ratio in RATIOS:
            if (n_cells, ratio) in done:
                print(f"skip complete cells={n_cells} ratio={ratio}")
                continue
            tmp = args.out.with_name(f"_ratio_tmp_{n_cells}_{ratio}.json")
            ns = types.SimpleNamespace(
                n_cells=n_cells, n_max=2, n_mc=n_mc, ntraj=ntraj,
                disorder_levels=[0.0], wigner_grid=31,
                initial_state_mode="edge_site", kerr=0.5,
                gamma_edge=ratio * GAMMA_BULK, gamma_bulk=GAMMA_BULK,
                t_stop=rcw.PERIOD, n_time=9, seed=20260814, output=tmp)
            print(f"cells={n_cells} ratio={ratio} ...", flush=True)
            rcw.run_dynamics(ns)
            recs = [r for r in json.loads(tmp.read_text())["records"]
                    if r["loss_profile"] == "edge_low_loss"]
            mean, se, n = paired_contrast(recs)
            matrix.append({"n_cells": n_cells, "ratio": ratio,
                           "contrast_mean": mean, "contrast_se": se, "n": n,
                           "sign": int((mean > 0) - (mean < 0))})
            args.out.write_text(json.dumps({"cells": matrix}, indent=2) + "\n")
            print(f"  contrast {mean:+.6f} +- {se:.2e} (n={n})")
    print(f"saved {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
