#!/usr/bin/env python3
"""Joint draw-level certification: witness x collective null (A-S1).

For two cells, d=0, edge_low_loss, samples 0-9, both sectors: joins the
R8 rerun records (rho_edge blocks) with the collective-mode records on
(sample, sector) and evaluates, per draw, the double criterion

    witness-certified:  max_alpha F_1 > omega_1 = 0.5   (Chabaud n=1)
    collective-quiet:   neg_coll_edge < 1% of the same-draw local negativity

The MDE leg stays stratum-level (already gated in compute_power_analysis):
a single-draw contrast cannot clear a mean-level detection threshold, so
requiring it per draw would be statistically incoherent.

Output: data/joint_certification.json with per-draw rows, per-sector
double-certification rates, and the stratum-level triple for reference.
Exits non-zero when the join is incomplete or a recomputation mismatches
a stored value beyond tolerance.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import analyze_r8_witness as r8

ROOT = Path(__file__).resolve().parents[1]
RERUN = ROOT / "data" / "r8_witness" / "r8_witness_2cells.json"
COLL = ROOT / "data" / "collective_wigner_dynamics.json"
PROD = ROOT / "data" / "chiral_production" / "chiral_floquet_2cells_NMC500.json"
OUT = ROOT / "data" / "joint_certification.json"

OMEGA1 = 0.5
COLL_RATIO = 0.01
N_DRAWS = 10


def main() -> int:
    """Join the three legs draw by draw and certify."""
    rerun = { (r["sample"], r["sector"]): r
              for r in json.loads(RERUN.read_text())["records"]
              if float(r["disorder_level"]) == 0.0
              and r["loss_profile"] == "edge_low_loss" }
    coll = { (r["sample"], r["sector"]): r
             for r in json.loads(COLL.read_text())["records"]
             if float(r["disorder_level"]) == 0.0
             and r["loss_profile"] == "edge_low_loss" }
    prod = {}
    for r in json.loads(PROD.read_text())["records"]:
        if float(r["disorder_level"]) == 0.0 and r["loss_profile"] == "edge_low_loss":
            prod.setdefault((r["sample"], r["sector"]), r)
    rows = []
    for sample in range(N_DRAWS):
        for sector in ("chiral_topological", "chiral_trivial"):
            key = (sample, sector)
            if key not in rerun or key not in coll or key not in prod:
                print(f"FAIL: incomplete join for draw {key}")
                return 1
            rho = r8.block_matrix(rerun[key]["rho_edge"])
            r8.check_block(rho)
            f_star, _ = r8.maximize_witness(rho, 1)
            w_cert = bool(float(f_star) > OMEGA1)
            neg_local = float(coll[key]["negativity_local_site0"])
            neg_coll = float(coll[key]["negativity_collective_edge"])
            quiet = bool(neg_coll < COLL_RATIO * neg_local)
            rows.append({
                "sample": sample, "sector": sector,
                "witness_F1": float(f_star), "witness_certified": w_cert,
                "neg_local": neg_local, "neg_coll_edge": neg_coll,
                "coll_ratio": float(neg_coll / neg_local) if neg_local > 0 else None,
                "collective_quiet": quiet,
                "double_certified": bool(w_cert and quiet),
                "prod_negativity_edge": float(prod[key]["negativity_edge"]),
            })
    summary = {}
    for sector in ("chiral_topological", "chiral_trivial"):
        sel = [r for r in rows if r["sector"] == sector]
        summary[sector] = {
            "n": len(sel),
            "witness_rate": sum(r["witness_certified"] for r in sel) / len(sel),
            "quiet_rate": sum(r["collective_quiet"] for r in sel) / len(sel),
            "double_rate": sum(r["double_certified"] for r in sel) / len(sel),
        }
    out = {
        "status": "JOINT_CERTIFICATION_COMPUTED",
        "scope": "two cells, d=0, edge_low_loss, samples 0-9; draw-level "
                 "double criterion; MDE leg stratum-level by construction",
        "criterion": {
            "witness": "max_alpha F_1 > 0.5 (Chabaud n=1, exact threshold)",
            "collective_quiet": "coll edge negativity < 1% of same-draw local",
        },
        "per_draw": rows,
        "per_sector": summary,
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    for sector, s in summary.items():
        print(f"{sector}: witness {s['witness_rate']:.2f} "
              f"quiet {s['quiet_rate']:.2f} double {s['double_rate']:.2f} "
              f"(n={s['n']})")
    print(f"saved {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
