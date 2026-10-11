#!/usr/bin/env python3
"""Fock-cutoff convergence check: N_max=4 vs production N_max=2 at four cells.

Audit follow-up item "4-cell N_max=4".  The production four-cell ensemble
(``data/chiral_production/chiral_floquet_4cells_NMC500.json``) truncates
each site at two Fock states inside the exact Q<=2 subspace; this script
re-simulates the first 100 disorder draws of every (profile, disorder)
stratum with ``--n-max 4`` (exact Q<=4 subspace, dim 495 of 5^8) under
identical seeds and computes paired per-record differences plus the
sector-paired Wigner-contrast shift.

Gates (script aborts on failure):
  1. all 8000 rerun records find a production partner;
  2. every stratum's sector-paired contrast shift has a 95% upper bound
     ``|mean| + 1.96*SE`` strictly below the operational floor
     FLOOR = 2e-3: at 95% confidence the N_max=4 truncation can never
     move any stratum contrast by a floor-scale amount.  Means, SEs and
     sigma-levels are reported so the reader sees which strata have
     resolved (but sub-floor) shifts.

Output:
    data/nmax4_cutoff.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FLOOR = 2e-3  # operational +-floor of the contrast criterion (dimensionless)


def main() -> None:
    """Pair the rerun with production and report cutoff-shift statistics."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path,
                   default=ROOT / "data" / "nmax4_cutoff.json")
    args = p.parse_args()

    prod = json.loads((ROOT / "data" / "chiral_production"
                       / "chiral_floquet_4cells_NMC500.json").read_text())
    new = json.loads((ROOT / "data" / "chiral_cutoff_nmax4"
                      / "chiral_floquet_4cells_NMC500_nmax4.json").read_text())
    pkey = {(r["loss_profile"], r["disorder_level"], r["sample"],
             r["sector"]): r for r in prod["records"]}

    pairs, missing = [], 0
    for r in new["records"]:
        k = (r["loss_profile"], r["disorder_level"], r["sample"], r["sector"])
        q = pkey.get(k)
        if q is None:
            missing += 1
            continue
        pairs.append((r, q))
    gates = {"all_paired": missing == 0 and len(pairs) == len(new["records"]),
             "shift_upper95_below_floor": True}

    per_record = {}
    for field in ("negativity_edge", "negativity_center", "edge_population"):
        d = np.array([a[field] - b[field] for a, b in pairs])
        per_record[field] = {"mean": float(d.mean()), "se": float(d.std(ddof=1) / np.sqrt(len(d))),
                             "max_abs": float(np.abs(d).max())}

    strata = {}
    for profile in ("edge_low_loss", "uniform_loss"):
        for disorder in sorted({a["disorder_level"] for a, _ in pairs}):
            sel = [(a, b) for a, b in pairs
                   if a["loss_profile"] == profile
                   and a["disorder_level"] == disorder]
            by_sector = {}
            for sector in ("chiral_topological", "chiral_trivial"):
                g = [(a, b) for a, b in sel if a["sector"] == sector]
                dn = np.array([a["negativity_edge"] - b["negativity_edge"]
                               for a, b in g])
                by_sector[sector] = {"mean_dn": float(dn.mean()),
                                     "se_dn": float(dn.std(ddof=1) / np.sqrt(len(dn)))}
            # sector-paired contrast: (topo - trivial), paired over samples
            tt = {s: {a["sample"]: (a, b) for a, b in sel if a["sector"] == s}
                  for s in ("chiral_topological", "chiral_trivial")}
            samples = sorted(set(tt["chiral_topological"]) & set(tt["chiral_trivial"]))
            c_new = np.array([tt["chiral_topological"][s][0]["negativity_edge"]
                              - tt["chiral_trivial"][s][0]["negativity_edge"]
                              for s in samples])
            c_old = np.array([tt["chiral_topological"][s][1]["negativity_edge"]
                              - tt["chiral_trivial"][s][1]["negativity_edge"]
                              for s in samples])
            shift = c_new - c_old
            se = float(shift.std(ddof=1) / np.sqrt(len(shift)))
            mean = float(shift.mean())
            strata[f"{profile}_d{disorder:g}"] = {
                "n": len(samples),
                "contrast_nmax4": float(c_new.mean()),
                "contrast_nmax2": float(c_old.mean()),
                "shift_mean": mean,
                "shift_se": se,
                "shift_sigma": mean / se if se > 0 else 0.0,
                "shift_upper95": abs(mean) + 1.96 * se,
            }
            if abs(mean) + 1.96 * se >= FLOOR:
                gates["shift_upper95_below_floor"] = False

    out = {
        "status": ("NMAX4_CUTOFF_COMPLETE" if all(gates.values())
                   else "NMAX4_CUTOFF_FAILED"),
        "claim_scope": "Fock-truncation convergence of the four-cell "
                       "production ensemble (N_max=4, exact Q<=4 subspace, "
                       "paired seeds, first 100 draws per stratum)",
        "not_supported": ["topological protection theorem",
                          "experimental validation", "metrological gain"],
        "design": {"n_cells": 4, "n_max_new": 4, "n_max_prod": 2,
                   "subspace_dim_new": 495, "subspace_dim_prod": 45,
                   "n_mc": 500, "ntraj": 10, "floor": FLOOR},
        "gates": gates,
        "per_record": per_record,
        "strata": strata,
        "software": {"numpy": np.__version__},
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(f"gates: {gates}")
    for f, v in per_record.items():
        print(f"  {f}: mean={v['mean']:+.3e} se={v['se']:.1e} max|d|={v['max_abs']:.3e}")
    for k, v in strata.items():
        print(f"  {k}: shift={v['shift_mean']:+.3e} +/- {v['shift_se']:.1e} "
              f"({v['shift_sigma']:+.1f} sigma, upper95={v['shift_upper95']:.1e}) "
              f"(c4={v['contrast_nmax4']:+.6f} c2={v['contrast_nmax2']:+.6f})")
    print(f"saved {args.output}")
    if not all(gates.values()):
        raise SystemExit("FAIL: " + ", ".join(k for k, v in gates.items()
                                              if not v))


if __name__ == "__main__":
    main()
