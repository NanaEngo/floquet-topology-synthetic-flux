#!/usr/bin/env python3
"""Plot the FL_QOM phase-scan audit.

The figure distinguishes framework winding output from the independent gap
Diagnostic. It is an audit figure, not a claim of a validated open-system
invariant phase diagram.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "floquet_phase_audit.json"
OUT_PDF = ROOT / "figures" / "floquet_phase_audit.pdf"
OUT_PNG = ROOT / "figures" / "floquet_phase_audit.png"


def main() -> None:
    data = json.loads(DATA.read_text())
    records = data["records"]
    j1 = sorted({r["J1"] for r in records})
    j2 = sorted({r["J2"] for r in records})
    shape = (len(j1), len(j2))
    index = {(r["J1"], r["J2"]): r for r in records}
    w20 = np.array([[index[(x, y)]["framework_winding_n20"] for y in j2] for x in j1])
    w40 = np.array([[index[(x, y)]["framework_winding_n40"] for y in j2] for x in j1])
    gap = np.array([[index[(x, y)]["min_quasienergy_gap_n40"] for y in j2] for x in j1])
    unstable = w20 != w40

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9), constrained_layout=True)
    extent = [j2[0], j2[-1], j1[0], j1[-1]]
    cmap = "viridis"
    for ax, values, letter, label in [
        (axes[0], w20, "(a)", "reported value"),
        (axes[1], w40, "(b)", "reported value"),
    ]:
        im = ax.imshow(values, origin="lower", extent=extent, aspect="auto", cmap=cmap)
        ax.set_xlabel(r"$J_2$")
        ax.set_ylabel(r"$J_1$")
        ax.text(0.02, 0.98, letter, transform=ax.transAxes,
                va="top", ha="left", fontsize=11)
        fig.colorbar(im, ax=ax, label=label)
        # Highlight the source-model point used in the N=6 rerun.
        ax.plot(1.0, 0.1, marker="*", markersize=9, color="white", markeredgecolor="black")
        ax.plot(1.0, 1.0, marker="s", markersize=6, color="red", markeredgecolor="black")

    im = axes[2].imshow(gap, origin="lower", extent=extent, aspect="auto", cmap="magma")
    axes[2].contour(*np.meshgrid(j2, j1), unstable.astype(float), levels=[0.5], colors="cyan", linewidths=1.2)
    axes[2].set_xlabel(r"$J_2$")
    axes[2].set_ylabel(r"$J_1$")
    axes[2].text(0.02, 0.98, "(c)", transform=axes[2].transAxes,
                 va="top", ha="left", fontsize=11)
    fig.colorbar(im, ax=axes[2], label="phase separation (rad)")
    axes[2].plot(1.0, 0.1, marker="*", markersize=9, color="white", markeredgecolor="black")
    axes[2].plot(1.0, 1.0, marker="s", markersize=6, color="red", markeredgecolor="black")

    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=300)
    print(f"Saved {OUT_PDF}")
    print(f"Saved {OUT_PNG}")
    print(f"resolution-disagreement points: {int(unstable.sum())}/{unstable.size}")


if __name__ == "__main__":
    main()
