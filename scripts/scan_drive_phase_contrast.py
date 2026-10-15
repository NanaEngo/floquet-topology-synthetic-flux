#!/usr/bin/env python3
"""Drive-phase contrast scan for the Track A sign pattern (A-S3).

Reuses run_collective_wigner.run_dynamics at fractional stroboscopic
times t_stop in {0.25, 0.5, 0.75, 1.0} x PERIOD (two cells, d=0, both
loss profiles, edge sector pair) and records the paired
topological-minus-trivial site-0 contrast per phase. Tests whether the
reported sign pattern is a stroboscopic-t=T accident. HPC-scale at full
sampling; --smoke runs one phase. Output: data/drive_phase_contrast.json.
"""
from __future__ import annotations

import argparse
import json
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_collective_wigner as rcw

OUT = ROOT / "data" / "drive_phase_contrast.json"
PHASES = (0.25, 0.5, 0.75, 1.0)


def main() -> int:
    """Run the drive-phase scan with per-phase checkpointing."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    phases = (1.0,) if args.smoke else PHASES
    n_mc, ntraj = (2, 4) if args.smoke else (10, 10)
    try:
        rows = json.loads(args.out.read_text())["phases"]
        done = {r["phase"] for r in rows}
    except (FileNotFoundError, KeyError):
        rows, done = [], set()
    for phase in phases:
        if phase in done:
            print(f"skip complete phase={phase}")
            continue
        tmp = args.out.with_name(f"_phase_tmp_{phase}.json")
        ns = types.SimpleNamespace(
            n_cells=2, n_max=2, n_mc=n_mc, ntraj=ntraj,
            disorder_levels=[0.0], wigner_grid=31,
            initial_state_mode="edge_site", kerr=0.5,
            gamma_edge=1e-4, gamma_bulk=0.2,
            t_stop=phase * rcw.PERIOD, n_time=9, seed=20260814,
            output=tmp)
        print(f"phase={phase} ...", flush=True)
        rcw.run_dynamics(ns)
        recs = json.loads(tmp.read_text())["records"]
        entry = {"phase": phase, "profiles": {}}
        import numpy as np
        for profile in ("edge_low_loss", "uniform_loss"):
            groups = {}
            for r in recs:
                if r["loss_profile"] != profile:
                    continue
                groups.setdefault((r["sample"], r["disorder_level"]), {})[r["sector"]] = r
            vals = np.array([g["chiral_topological"]["negativity_local_site0"]
                             - g["chiral_trivial"]["negativity_local_site0"]
                             for g in groups.values() if len(g) == 2])
            entry["profiles"][profile] = {
                "mean": float(vals.mean()), "se": float(vals.std(ddof=1) / len(vals) ** 0.5),
                "n": len(vals)} if len(vals) else {"mean": 0.0, "se": 0.0, "n": 0}
        rows.append(entry)
        args.out.write_text(json.dumps({"phases": rows}, indent=2) + "\n")
        print(f"  edge {entry['profiles']['edge_low_loss']['mean']:+.6f} "
              f"uniform {entry['profiles']['uniform_loss']['mean']:+.6f}")
    print(f"saved {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
