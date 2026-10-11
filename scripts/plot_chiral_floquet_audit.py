#!/usr/bin/env python3
"""Plot the chiral Floquet phase diagram and finite-chain edge states."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
plt.style.use(["seaborn-v0_8-ticks", "seaborn-v0_8-colorblind"])
import numpy as np

from matplotlib import patheffects as path_effects

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "chiral_floquet_phase_audit.json"
OUT_PDF = ROOT / "figures" / "chiral_floquet_phase_audit.pdf"
OUT_PNG = ROOT / "figures" / "chiral_floquet_phase_audit.png"


def matrix_for(records, key):
    xs = sorted({r["ja_intra_step_A"] for r in records})
    ys = sorted({r["jb_inter_step_B"] for r in records})
    idx = {(r["ja_intra_step_A"], r["jb_inter_step_B"]): r for r in records}
    z = np.array([[idx[(x, y)]["resolution_results"][-1][key] for y in ys] for x in xs])
    return xs, ys, z


def main():
    d = json.loads(DATA.read_text())
    records = d["scan"]["records"]
    xs, ys, nu0 = matrix_for(records, "nu_0")
    _, _, nupi = matrix_for(records, "nu_pi")
    edge = np.array([[r["edge_diagnostics"]["12"]["n_edge_near_zero"] for r in records
                      if r["ja_intra_step_A"] == x and r["jb_inter_step_B"] == y][0]
                     for x in xs for y in ys]).reshape(len(xs), len(ys))
    extent = [ys[0], ys[-1], xs[0], xs[-1]]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    for ax, values, letter, cmap in [
        (axes[0], nu0, "(a)", "coolwarm"),
        (axes[1], nupi, "(b)", "RdBu"),
        (axes[2], edge, "(c)", "viridis"),
    ]:
        im = ax.imshow(values, origin="lower", extent=extent, aspect="auto", cmap=cmap)
        ax.set_xlabel(r"$J_B$ (intercell pulse)")
        ax.set_ylabel(r"$J_A$ (intracell pulse)")
        ax.text(0.98, 0.98, letter, transform=ax.transAxes,
                va="top", ha="right", fontsize=11, color="white",
                path_effects=[path_effects.withStroke(linewidth=1.5, foreground="black")])
        fig.colorbar(im, ax=ax)
        ax.plot(1.2, 0.1, marker="*", markersize=9, color="white", markeredgecolor="black",
                label="topological point")
        ax.plot(0.1, 1.2, marker="s", markersize=6, color="cyan", markeredgecolor="black",
                label="trivial point")
    axes[0].legend(loc="lower left", frameon=True, fontsize=8)
    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=300)
    print(f"Saved {OUT_PDF}")
    print(f"Saved {OUT_PNG}")


if __name__ == "__main__":
    main()
