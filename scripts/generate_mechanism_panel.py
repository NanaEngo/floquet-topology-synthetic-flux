#!/usr/bin/env python3
"""Generate the main-text mechanism panel (Fig. 4) for Track A.

Three panels:
  (a) topological edge-mode profiles |u_{e,j}|^2 at 2, 3, 4 cells;
  (b) measured decay rate Gamma_meas vs the mode-weighted rate Gamma_eff;
      the two-photon dwell is shown only as a size-specific diagnostic;
  (c) per-size kappa_N = gamma_eb/Lambda_2^coh for topological and trivial
      sectors, compared with the exact single-site reference kappa_single=3.411.

Regenerated 2026-08-15 after the independent Lambda_2 recomputation: the
previous universal-kappa panel was removed because the archived Lambda_2 values
were not reproducible from the documented definition. Since 2026-10-02 the
dwell integrals and the trivial-sector prefactor are read from the
self-checking derivation archive (data/closedform_kappa_derivation.json,
produced by scripts/derive_closedform_kappa.py) instead of hardcoded literals.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.style.use(["seaborn-v0_8-ticks", "seaborn-v0_8-colorblind"])
import numpy as np

from run_chiral_floquet_dynamics import (
    CANDIDATES,
    edge_mode_vector,
    rate_vector,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "mechanism_panel.pdf"

N_CELLS_LIST = (2, 3, 4)
KAPPA_SINGLE = 3.411      # exact bare-site reference, SI Eq. S20
KAPPA_ASYMP = 2.0         # 2|alpha|^2 with alpha=1

# Rate-weighted two-photon dwell Lambda_2^coh (1/period), edge-low-loss, zero
# disorder — read from the self-checking derivation archive so the panel
# cannot drift from the documented definition.
_DERIVED = json.loads(
    (ROOT / "data" / "closedform_kappa_derivation.json").read_text()
)["sizes"]
LAMBDA2_COH = {int(n): s["lambda2_coh"] for n, s in _DERIVED.items()}

# Production-parameter rates read from the participation-bound artifact so the
# panel cannot drift from the archived fits. Tuple: Gamma_eff, Gamma_meas,
# Lambda_2^coh, gamma_eb = Gamma_meas - Gamma_eff, kappa_N = gamma_eb/Lambda_2^coh.
# All rates are in units of 1/period (period = 2 pi).
_PARTICIPATION = ROOT / "data" / "participation_bound_verification.json"


def load_meas():
    """Load per-size production rates from `participation_bound_verification.json`.

    Returns:
        dict[int, tuple[float, float, float, float, float]]: keys are cell
        counts 2-4; values are (Gamma_eff, Gamma_meas, Lambda_2^coh, gamma_eb,
        kappa_N) in 1/period, except the dimensionless kappa_N.
    """
    bound = json.loads(_PARTICIPATION.read_text())["bound_verification"]
    meas = {}
    for row in bound:
        n = row["n_cells"]
        geb = row["gamma_meas"] - row["gamma_eff"]
        l2 = LAMBDA2_COH[n]
        meas[n] = (row["gamma_eff"], row["gamma_meas"], l2, geb, geb / l2)
    return meas


MEAS = load_meas()
# Trivial-sector effective prefactor kappa_eff = gamma_eb/Lambda_2 from the
# same derivation archive (withdrawn 1.88/1.78/1.72 literals replaced 2026-10-02).
TRIVIAL_KAPPA = {int(n): s["trivial"]["kappa_eff"] for n, s in _DERIVED.items()}

fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.1))
colors = ["#1f77b4", "#d62728", "#2ca02c"]

# ---- (a) edge-mode profiles ----
ax = axes[0]
for n, c in zip(N_CELLS_LIST, colors):
    u = edge_mode_vector(n, CANDIDATES["chiral_topological"])
    w = np.abs(u) ** 2
    ax.plot(np.arange(2 * n), w, "o-", ms=4, lw=1.2, color=c, label=f"{n} cells")
ax.set_xlabel(r"site $j$")
ax.set_ylabel(r"$|u_{e,j}|^2$")
ax.set_ylim(0, 0.55)
ax.legend(fontsize=7, frameon=False)
ax.text(0.02, 0.98, "(a)", transform=ax.transAxes,
        va="top", ha="left", fontsize=11)

# ---- (b) decay rates ----
ax = axes[1]
ns = np.arange(len(N_CELLS_LIST))
for i, n in enumerate(N_CELLS_LIST):
    ge, gm, l2, geb, kn = MEAS[n]
    ax.plot([i - 0.12, i], [ge, gm], "o-", ms=5, lw=0.9, color=colors[i])
    ax.text(i, gm + 0.05, rf"{gm:.3f}", ha="center", fontsize=7, color=colors[i])
    ax.text(i - 0.12, ge + 0.05, rf"{ge:.3f}", ha="center", fontsize=6, color=colors[i])
ax.set_xticks(ns)
ax.set_xticklabels([f"{n} cells" for n in N_CELLS_LIST], fontsize=8)
ax.set_ylabel(r"decay rate $\Gamma$")
ax.text(0.02, 0.98, "(b)", transform=ax.transAxes,
        va="top", ha="left", fontsize=11)
ax.set_ylim(0, 0.85)
ax.annotate(r"$\Gamma_{\rm meas}$", xy=(0, 0.33), xytext=(0.25, 0.40),
            fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.7))
ax.annotate(r"$\Gamma_{\rm eff}$", xy=(-0.12, 0.185), xytext=(0.18, 0.235),
            fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.7))

# ---- (c) size-specific efficiency ----
ax = axes[2]
positions = np.arange(len(N_CELLS_LIST))
topo_k = np.array([MEAS[n][4] for n in N_CELLS_LIST])
triv_k = np.array([TRIVIAL_KAPPA[n] for n in N_CELLS_LIST])
ax.plot(positions - 0.08, topo_k, "o-", ms=5, lw=0.9, color=colors[0], label=r"topological $\kappa_N$")
ax.plot(positions + 0.08, triv_k, "^--", ms=5, lw=0.9, color="0.4", label=r"trivial $\kappa_N$")
ax.axhline(KAPPA_SINGLE, ls="-", lw=1.0, color="k", label=r"$\kappa_{\rm single}=3.411$")
ax.axhline(KAPPA_ASYMP, ls=":", lw=1.0, color="0.6", label=r"$2|\alpha|^2=2.0$")
ax.set_xlabel("chain size")
ax.set_xticks(positions)
ax.set_xticklabels([f"{n}" for n in N_CELLS_LIST])
ax.set_ylabel(r"per-size $\kappa_N$")
ax.set_ylim(0, 13)
ax.legend(fontsize=6.2, frameon=False, loc="upper right")
ax.text(0.02, 0.98, "(c)", transform=ax.transAxes,
        va="top", ha="left", fontsize=11)

for ax in axes:
    ax.tick_params(labelsize=8)

fig.tight_layout()
fig.savefig(OUT, bbox_inches="tight")
fig.savefig(ROOT / "figures" / "mechanism_panel.png", dpi=200, bbox_inches="tight")
print(f"saved {OUT}")
