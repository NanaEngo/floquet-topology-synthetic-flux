#!/usr/bin/env python3
"""Plot the micromotion-covariant edge-mode dynamics."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
plt.style.use(["seaborn-v0_8-ticks", "seaborn-v0_8-colorblind"])
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "edge_mode_covariant_analysis.json"
OUT = ROOT / "figures" / "edge_mode_covariant_summary.pdf"
OUT_PNG = ROOT / "figures" / "edge_mode_covariant_summary.png"


def main() -> None:
    data = json.loads(INPUT.read_text())
    lookup = {}
    for item in data["files"]:
        name = Path(item["input"]).name
        regime = "Kerr" if "kerr" in name else "passive"
        prep = "common site" if "common" in name else "edge mode"
        for summary in item["summaries"]:
            if summary["loss_profile"] == "edge_low_loss" and summary["disorder_level"] == 0.0:
                m = summary["edge_mode_population_topological_minus_trivial"]
                lookup[(prep, regime)] = (m["mean"], m["ci95_low"], m["ci95_high"])

    labels = [("common site", "passive"), ("common site", "Kerr"),
              ("edge mode", "passive"), ("edge mode", "Kerr")]
    means = np.array([lookup[k][0] for k in labels])
    lower = np.array([lookup[k][0] - lookup[k][1] for k in labels])
    upper = np.array([lookup[k][2] - lookup[k][0] for k in labels])

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    x = np.arange(len(labels))
    colors = ["#2166ac", "#67a9cf", "#b2182b", "#ef8a62"]
    ax.errorbar(x, means, yerr=np.vstack([lower, upper]), fmt="none", ecolor="#222222",
                elinewidth=1.2, capsize=3, zorder=3)
    ax.bar(x, means, color=colors, edgecolor="#222222", linewidth=0.5, width=0.68, zorder=2)
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_xticks(x, ["common\npassive", "common\nKerr", "edge-mode\npassive", "edge-mode\nKerr"])
    ax.tick_params(axis="x", labelsize=8)
    ax.set_ylabel(r"$\Delta n_e$ (topological $-$ trivial)")
    ax.text(0.02, 0.98, "(a)", transform=ax.transAxes,
            va="top", ha="left", fontsize=11)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6, zorder=0)

    ax = axes[1]
    residuals = [0.039564269521941545, 0.38257089148882073]
    ax.bar([0, 1], residuals, color=["#2166ac", "#b2182b"], edgecolor="#222222", linewidth=0.5)
    ax.set_xticks([0, 1], ["common\nsite", "edge-mode\nprep."])
    ax.set_ylabel(r"max $|R_e(t)|$")
    ax.text(0.02, 0.98, "(b)", transform=ax.transAxes,
            va="top", ha="left", fontsize=11)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6, zorder=0)
    ax.set_ylim(0, max(residuals) * 1.18)

    fig.tight_layout(pad=1.0)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(OUT)
    print(OUT_PNG)


if __name__ == "__main__":
    main()
