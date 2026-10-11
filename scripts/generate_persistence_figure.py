#!/usr/bin/env python3
"""Rebuild Track A's persistence display from the archived trace arrays.

The arrays are treated as numerical evidence, not as proof of a theorem.  The script
records their provenance, time grid, sample sizes, means, and pointwise 95% intervals.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT.parents[1] / "code" / "_cache"
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
DATA.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

FILES = {
    "W_plus_one_half": CACHE / "persistence_mc_W+1_2.npy",
    "W_minus_one_half": CACHE / "persistence_mc_W-1_2.npy",
    "W_zero": CACHE / "persistence_mc_W0.npy",
}

traces: dict[str, np.ndarray] = {}
for label, path in FILES.items():
    if not path.exists():
        raise FileNotFoundError(path)
    arr = np.load(path)
    if arr.ndim != 2 or arr.shape[0] != 500:
        raise ValueError(f"Unexpected shape for {path}: {arr.shape}")
    traces[label] = np.asarray(arr, dtype=float)

n_time = next(iter(traces.values())).shape[1]
time = np.linspace(0.0, 10.0, n_time)
summary = {
    "source_files": {k: str(v.relative_to(ROOT.parents[1])) for k, v in FILES.items()},
    "source_kind": "archived numerical trace arrays; no new quantum propagation was run",
    "n_realizations": 500,
    "n_time_points": int(n_time),
    "time_grid": {"start": 0.0, "stop": 10.0, "n": int(n_time)},
    "interpretation": "descriptive finite-time persistence comparison; not a theorem or asymptotic bound",
    "groups": {},
}

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(6.2, 3.8))
colors = {"W_plus_one_half": "#2166ac", "W_minus_one_half": "#d6604d", "W_zero": "#4d4d4d"}
labels = {
    "W_plus_one_half": r"archived label $W=+1/2$",
    "W_minus_one_half": r"archived label $W=-1/2$",
    "W_zero": r"archived label $W=0$",
}

for key, arr in traces.items():
    mean = arr.mean(axis=0)
    lo, hi = np.percentile(arr, [2.5, 97.5], axis=0)
    summary["groups"][key] = {
        "mean_initial": float(mean[0]),
        "mean_final": float(mean[-1]),
        "mean_at_nearest_t_1": float(mean[np.argmin(abs(time - 1.0))]),
        "nearest_time_to_1": float(time[np.argmin(abs(time - 1.0))]),
        "pointwise_95_percent_interval_at_nearest_t_1": [
            float(lo[np.argmin(abs(time - 1.0))]),
            float(hi[np.argmin(abs(time - 1.0))]),
        ],
    }
    ax.plot(time, mean, lw=2, color=colors[key], label=labels[key])
    ax.fill_between(time, lo, hi, color=colors[key], alpha=0.16, linewidth=0)

ax.set_xlabel(r"normalized time $t$")
ax.set_ylabel(r"integrated negative Wigner volume $\delta_{\rm neg}$")
ax.legend(frameon=False, fontsize=8)
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
fig.savefig(FIGURES / "trackA_persistence_rebuilt.pdf", bbox_inches="tight")
fig.savefig(FIGURES / "trackA_persistence_rebuilt.png", dpi=300, bbox_inches="tight")
summary["figure"] = {
    "pdf": "figures/trackA_persistence_rebuilt.pdf",
    "png": "figures/trackA_persistence_rebuilt.png",
}
(DATA / "trackA_persistence_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
