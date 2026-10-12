#!/usr/bin/env python3
"""Replay gate for the R8 witness chain (A-2).

Recomputes the n=1 Chabaud witness on a 3-draw subset (samples 0-2,
both sectors, two cells, d=0, edge_low_loss) through the same pure
functions analyze_r8_witness.py uses (block_matrix, check_block,
maximize_witness) and checks each sector subset mean against the
published r8_witness_summary.json stratum. Tolerance scales with the
subset size (3 se*sqrt(N/n)) since a 3-draw mean carries its own
sampling noise; a systematic recomputation error would exceed it.

This exercises the analysis code path end to end (load -> parity-shaped
records -> witness maximization -> aggregation) without rerunning the
full record set. Exits non-zero on any mismatch.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import analyze_r8_witness as r8

ROOT = Path(__file__).resolve().parents[1]
RERUN = ROOT / "data" / "r8_witness" / "r8_witness_2cells.json"
SUMMARY = ROOT / "data" / "r8_witness_summary.json"
STRATUM = "2cells_d0_edge_low_loss"
N_REPLAY = 3
Z_TOL = 3.0

def main() -> int:
    """Replay the witness computation on a 3-draw subset."""
    rerun = json.loads(RERUN.read_text())
    summary = json.loads(SUMMARY.read_text())
    stratum = summary["strata"][STRATUM]
    n_pub = int(stratum["n"])
    fails = 0
    for sector in ("chiral_topological", "chiral_trivial"):
        sel = sorted(
            (r for r in rerun["records"]
             if r["sector"] == sector
             and float(r["disorder_level"]) == 0.0
             and r["loss_profile"] == "edge_low_loss"),
            key=lambda r: int(r["sample"]))[:N_REPLAY]
        if len(sel) < N_REPLAY:
            print(f"FAIL {sector}: only {len(sel)} records for replay")
            return 1
        vals = []
        for r in sel:
            rho = r8.block_matrix(r["rho_edge"])
            r8.check_block(rho)
            f_star, _ = r8.maximize_witness(rho, 1)
            vals.append(float(f_star))
        pub = summary["strata"][STRATUM]["witness"]["n1"][sector]
        mean = float(np.mean(vals))
        dev = abs(mean - pub["mean"])
        # subset-mean sampling noise: se*sqrt(N_pub/n_subset)
        tol = Z_TOL * pub["se"] * (n_pub / len(vals)) ** 0.5 + 1e-12
        ok = dev <= tol
        print(f"{'OK  ' if ok else 'FAIL'} {sector}: replay mean {mean:.4f} "
              f"vs published {pub['mean']:.4f} (dev {dev:.2e}, tol {tol:.2e}, "
              f"n={len(vals)} of {n_pub})")
        if not ok:
            fails += 1
    if fails:
        print(f"REPLAY FAILS: {fails}")
        return 1
    print("REPLAY FAILS: 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
