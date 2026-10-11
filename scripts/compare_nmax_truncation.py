"""Paired N_max=2 vs N_max=3 contrast comparison (same seeds, matched design).

For each system size (2, 3, 4 cells) and each (loss_profile, disorder_level),
the matched per-unit contrast delta = neg_edge(topological) - neg_edge(trivial)
is computed at both Fock cutoffs for the identical disorder realizations
(draw_seed and trajectory_seed must agree — asserted). Reports mean contrasts
per cutoff, the paired shift mean(delta_3 - delta_2) with bootstrap 95% CIs
(1e4 resamples, rng = default_rng(42)), per-unit sign agreement, and the
cross-cutoff correlation.

Output: data/nmax2_nmax3_paired_comparison.json

Usage: mamba run -n qom python scripts/compare_nmax_truncation.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROD2 = ROOT / "data" / "chiral_production"
PROD3 = ROOT / "data" / "chiral_production_nmax3"
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


def contrast_by_sample(path: Path, profile: str, disorder: float) -> dict[int, dict]:
    """Paired per-sample contrast keyed by sample, with seed provenance.

    Args:
        path: Production JSON file with a records list.
        profile (str): Loss profile to keep.
        disorder (float): Disorder level to keep.

    Returns:
        dict[int, dict]: sample -> {"delta", "seed", "topo", "triv"}.
    """
    data = json.loads(path.read_text())
    groups: dict[int, dict] = {}
    for r in data["records"]:
        if r["loss_profile"] == profile and float(r["disorder_level"]) == disorder:
            groups.setdefault(int(r["sample"]), {})[r["sector"]] = r
    out = {}
    for sample, pair in groups.items():
        if "chiral_topological" not in pair or "chiral_trivial" not in pair:
            continue
        topo, triv = pair["chiral_topological"], pair["chiral_trivial"]
        assert topo["draw_seed"] == triv["draw_seed"], "sector seeds differ"
        out[sample] = {
            "delta": topo["negativity_edge"] - triv["negativity_edge"],
            "seed": topo["draw_seed"],
            "topo": topo["negativity_edge"],
            "triv": triv["negativity_edge"],
        }
    return out


def bootstrap_mean_ci(vals: np.ndarray) -> list[float]:
    """Bootstrap-percentile 95% CI of the mean.

    Args:
        vals (np.ndarray): Independent unit values (n >= 2).

    Returns:
        list[float]: [2.5%, 97.5%] percentiles of the bootstrap mean.
    """
    resampled = np.array([
        vals[RNG.integers(0, vals.size, vals.size)].mean()
        for _ in range(N_RESAMPLE)
    ])
    return [float(x) for x in np.percentile(resampled, [2.5, 97.5])]


def best_file(prod: Path, cells: int) -> Path:
    cands = list(prod.glob(f"chiral_floquet_{cells}cells_NMC*.json"))
    return sorted(cands, key=lambda p: int(json.loads(p.read_text())["design"]["n_mc"]))[-1]


def main() -> None:
    out: dict = {
        "status": "PAIRED_TRUNCATION_COMPARISON_COMPUTED",
        "scope": ("matched N_max=2 vs N_max=3 per-unit contrast shift on identical "
                  "disorder realizations (seed-asserted); descriptive, bootstrap 95% CIs"),
        "design": {"n_resample": N_RESAMPLE, "rng": "default_rng(42)",
                   "observable": "negativity_edge(topological) - negativity_edge(trivial)",
                   "independent_unit": "disorder realization (paired across cutoffs)"},
        "sizes": {},
    }
    for cells in (2, 3, 4):
        path2, path3 = best_file(PROD2, cells), best_file(PROD3, cells)
        size_entry = {"nmax2_file": str(path2.relative_to(ROOT)),
                      "nmax3_file": str(path3.relative_to(ROOT)), "conditions": {}}
        for profile, disorder in CONDITIONS:
            d2 = contrast_by_sample(path2, profile, disorder)
            d3 = contrast_by_sample(path3, profile, disorder)
            common = sorted(set(d2) & set(d3))
            assert common, f"no common samples for {cells} cells {profile} d={disorder}"
            for s in common:
                assert d2[s]["seed"] == d3[s]["seed"], f"seed mismatch at sample {s}"
            v2 = np.array([d2[s]["delta"] for s in common])
            v3 = np.array([d3[s]["delta"] for s in common])
            shift = v3 - v2
            agree = np.mean(np.sign(v3) == np.sign(v2))
            r = (float(np.corrcoef(v2, v3)[0, 1])
                 if np.std(v2) > 0 and np.std(v3) > 0 else None)
            size_entry["conditions"][f"{profile}/d={disorder}"] = {
                "n": len(common),
                "mean_nmax2": float(v2.mean()),
                "ci95_nmax2": bootstrap_mean_ci(v2),
                "mean_nmax3": float(v3.mean()),
                "ci95_nmax3": bootstrap_mean_ci(v3),
                "paired_shift_mean": float(shift.mean()),
                "paired_shift_ci95": bootstrap_mean_ci(shift),
                "sign_agreement": float(agree),
                "pearson_r": r,
            }
        out["sizes"][f"{cells}cells"] = size_entry
    dest = ROOT / "data" / "nmax2_nmax3_paired_comparison.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dest}")
    for ckey, entry in out["sizes"].items():
        c = entry["conditions"]["edge_low_loss/d=0.0"]
        print(f"  {ckey} edge_low/d=0: shift={c['paired_shift_mean']:+.6f} "
              f"CI{[round(x, 6) for x in c['paired_shift_ci95']]} "
              f"agree={c['sign_agreement']:.3f} r={c['pearson_r'] if c['pearson_r'] is None else round(c['pearson_r'], 3)}")


if __name__ == "__main__":
    main()
