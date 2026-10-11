"""Regenerate data/bootstrap_contrast_cis.json from production JSONs.

Method: for each (cells, loss_profile, disorder_level), bootstrap-percentile 95% CI
over the independent disorder realizations of the matched chiral-minus-trivial
edge-negativity contrast, 1e4 resamples, rng = default_rng(42).

Usage:
    mamba run -n qom python scripts/compute_bootstrap_cis.py
    mamba run -n qom python scripts/compute_bootstrap_cis.py \
        --prod-dir data/chiral_production_nmax3 --out data/bootstrap_contrast_cis_nmax3.json
"""
import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(42)
N_RESAMPLE = 10_000
CONDITIONS = [
    ("edge_low_loss", 0.0),
    ("edge_low_loss", 0.5),
    ("edge_low_loss", 1.0),
    ("edge_low_loss", 1.5),
    ("uniform_loss", 0.0),
    ("uniform_loss", 0.5),
    ("uniform_loss", 1.0),
    ("uniform_loss", 1.5),
]


def paired_contrasts(path: Path, profile: str, disorder: float) -> list[float]:
    data = json.loads(path.read_text())
    groups = {}
    for r in data["records"]:
        if r["loss_profile"] == profile and r["disorder_level"] == disorder:
            groups.setdefault(int(r["sample"]), {})[r["sector"]] = r
    out = []
    for pair in groups.values():
        if "chiral_topological" in pair and "chiral_trivial" in pair:
            out.append(pair["chiral_topological"]["negativity_edge"]
                       - pair["chiral_trivial"]["negativity_edge"])
    return out


def best_file(prod: Path, cells: int) -> Path:
    cands = [p for p in prod.glob(f"chiral_floquet_{cells}cells_NMC*.json")]
    return sorted(cands, key=lambda p: int(json.loads(p.read_text())["design"]["n_mc"]))[-1]


def main() -> None:
    p = argparse.ArgumentParser(description="Bootstrap 95% CIs of the matched contrast")
    p.add_argument("--prod-dir", default="data/chiral_production",
                   help="production directory relative to ROOT")
    p.add_argument("--out", default="data/bootstrap_contrast_cis.json",
                   help="output JSON relative to ROOT")
    args = p.parse_args()
    prod = ROOT / args.prod_dir
    out = {}
    for cells in (2, 3, 4):
        path = best_file(prod, cells)
        out[f"{cells}cells"] = {}
        for profile, disorder in CONDITIONS:
            vals = np.asarray(paired_contrasts(path, profile, disorder), dtype=float)
            if vals.size == 0:
                raise SystemExit(f"no paired units for {cells} cells {profile} d={disorder}")
            resampled = np.array([
                vals[RNG.integers(0, vals.size, vals.size)].mean()
                for _ in range(N_RESAMPLE)
            ])
            out[f"{cells}cells"][f"{profile}/d={disorder}"] = {
                "n": int(vals.size),
                "mean": float(vals.mean()),
                "boot_ci95": [float(x) for x in np.percentile(resampled, [2.5, 97.5])],
            }
    dest = ROOT / args.out
    dest.write_text(json.dumps(out, indent=2) + "\n")
    ns = [out[f"{c}cells"]["edge_low_loss/d=0.0"]["n"] for c in (2, 3, 4)]
    print(f"wrote {dest} (cells 2/3/4, n_units {'/'.join(map(str, ns))})")


if __name__ == "__main__":
    main()