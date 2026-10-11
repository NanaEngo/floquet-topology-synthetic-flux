#!/usr/bin/env python3
"""Plot the completed deterministic N=6 finite-time source-model check."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
data = json.loads((ROOT / "data" / "trackA_n6_independent.json").read_text())
records = data["disorder_records"]
disorder = np.array([r["disorder_level"] for r in records])
topo = np.array([r["topological_mean"] for r in records])
topo_sd = np.array([r["topological_std"] for r in records])
reference = np.array([r["trivial_mean"] for r in records])
reference_sd = np.array([r["trivial_std"] for r in records])

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(6.2, 3.8))
ax.errorbar(disorder, topo, yerr=topo_sd, fmt="o-", lw=1.8, capsize=3,
            color="#2166ac", label=r"source case $(J_1,J_2)=(0.1,1.0)$")
ax.errorbar(disorder, reference, yerr=reference_sd, fmt="s--", lw=1.6, capsize=3,
            color="#4d4d4d", label=r"reference case $(J_1,J_2)=(1.0,1.0)$")
ax.set_xlabel("disorder amplitude")
ax.set_ylabel(r"endpoint $\delta_{\rm neg}(t=1)$")
ax.grid(axis="y", alpha=0.25)
ax.legend(frameon=False, fontsize=8)
ax.text(0.02, 0.98, "500 disorder samples; 10 trajectories/sample",
        transform=ax.transAxes, va="top", fontsize=8)
fig.tight_layout()
for suffix, kwargs in [("pdf", {}), ("png", {"dpi": 300})]:
    fig.savefig(ROOT / "figures" / f"trackA_n6_rebuild_summary.{suffix}",
                bbox_inches="tight", **kwargs)
print("saved trackA_n6_rebuild_summary.pdf/.png")
