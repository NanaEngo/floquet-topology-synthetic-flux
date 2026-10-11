#!/usr/bin/env python3
"""Robustness matrix (roadmap G11, Figure 5).

Panel (a): edge-low-loss contrast vs disorder, one curve per system size.
Panel (b): uniform-loss contrast vs disorder, same layout. The sign flip of the
edge-low-loss contrast between 2 and 3 cells, and the near-zero uniform-loss
contrast, are the prespecified negative findings.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
plt.style.use(["seaborn-v0_8-ticks", "seaborn-v0_8-colorblind"])
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "chiral_production"
OUT_PDF = ROOT / "figures" / "robustness_matrix.pdf"
OUT_PNG = ROOT / "figures" / "robustness_matrix.png"

FILES = {
    2: DATA / "chiral_floquet_2cells_NMC500.json",
    3: DATA / "chiral_floquet_3cells_NMC500.json",
    4: DATA / "chiral_floquet_4cells_NMC2000.json",
}


def load(cells):
    d = json.loads(FILES[cells].read_text())
    return d["records"]


def contrast_by(records, profile, disorder):
    """Matched sector contrast: mean(chiral) - mean(trivial) at one condition."""
    vals = {}
    for r in records:
        if r["loss_profile"] != profile or r["disorder_level"] != disorder:
            continue
        vals.setdefault(r["sector"], []).append(r["negativity_edge"])
    chiral = np.mean(vals["chiral_topological"])
    trivial = np.mean(vals["chiral_trivial"])
    # paired standard error of the contrast
    pairs = list(zip(vals["chiral_topological"], vals["chiral_trivial"]))
    diffs = np.asarray([c - t for c, t in pairs])
    return chiral - trivial, diffs.std(ddof=1) / np.sqrt(len(diffs)) if len(diffs) > 1 else 0.0


def main():
    records = {c: load(c) for c in (2, 3, 4)}
    disorders = [0.0, 0.5, 1.0, 1.5]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    colors = {2: "#1f77b4", 3: "#d62728", 4: "#2ca02c"}
    for ax, profile, letter in [
        (axes[0], "edge_low_loss", "(a)"),
        (axes[1], "uniform_loss", "(b)"),
    ]:
        ax.text(0.02, 0.98, letter, transform=ax.transAxes,
                va="top", ha="left", fontsize=11)
        for cells in (2, 3, 4):
            mu = []; se = []
            for d in disorders:
                c, e = contrast_by(records[cells], profile, d)
                mu.append(c); se.append(e)
            ax.errorbar(disorders, mu, yerr=se, marker="o", capsize=3, lw=1.6,
                        color=colors[cells], label=f"{cells} cells / {2 * cells} modes")
        ax.axhline(0.0, color="k", lw=0.8, ls=":")
        ax.set_xlabel("Disorder level")
        ax.set_ylabel(r"Contrast $\Delta$ (chiral $-$ trivial)")
        ax.legend(fontsize=8, frameon=False)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=300)
    print(f"Saved {OUT_PDF}")
    print(f"Saved {OUT_PNG}")


if __name__ == "__main__":
    main()
