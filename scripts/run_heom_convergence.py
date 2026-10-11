#!/usr/bin/env python3
"""HEOM hierarchy-depth convergence on the production Floquet-Lindblad model.

The archived HEOM record (``data/trackA_heom_validation.json``) is a
depth-1 check on a different (static, three-site) setup with no producer
script.  This script is the production-model counterpart: it integrates the
two-cell driven Kerr chain with every site coupled to a Drude-Lorentz bath
(``qutip.solver.heom``) at hierarchy depths 1/2/3, plus a Matsubara-terms
sensitivity point (Nk=3 at depth 1), and compares the terminal edge-site
Wigner negativity across depths.

Scope (deliberate): this is a convergence check of the HEOM solver
machinery on the production Hamiltonian, not a re-validation of the
production Lindblad runs -- the microscopic Drude-Lorentz bath differs
from the Markovian collapse operators, so an absolute offset against
``mesolve`` is expected and reported as context, not gated.

Gates (script aborts on failure):
  1. every (sector, depth) solve completes with a finite physical state;
  2. depth increments shrink: |neg(d3) - neg(d2)| <= |neg(d2) - neg(d1)|
     (hierarchy convergence; deterministic solver, no sampling noise).

Output:
    data/heom_convergence.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import qutip as qt
from qutip.solver.heom import DrudeLorentzBath, HEOMSolver

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_chiral_floquet_dynamics import (  # noqa: E402
    CANDIDATES,
    PERIOD,
    build_ops,
    driven_hamiltonian,
    initial_state,
    negativity,
)

ROOT = Path(__file__).resolve().parents[1]

# Bath operating point (dimensionless, drive-frequency units): strong,
# explicitly non-Markovian-corrected damping so that hierarchy depth
# matters.  lam/gamma_cutoff chosen so the Markov-equivalent rate is
# order-unity (far above the production edge rate); T/Nk mirror the
# archived check (T=0.5, Nk=2) with one Nk=3 sensitivity point.
LAM = 0.3
GAMMA_CUTOFF = 1.0
T_BATH = 0.5
NK = 2
NK_SENSITIVITY = 3
DEPTHS = (1, 2, 3)


def run_point(n_cells: int, sector: str, depth: int, nk: int) -> dict:
    """One HEOM solve: terminal edge/centre negativity and populations.

    Args:
        n_cells (int): Unit cells (only two supported: density-matrix scale).
        sector (str): ``chiral_topological`` or ``chiral_trivial``.
        depth (int): HEOM hierarchy maximum depth (dimensionless).
        nk (int): Matsubara expansion terms per bath (dimensionless).

    Returns:
        dict: ``negativity_edge``, ``negativity_center`` (dimensionless),
        ``edge_population`` (dimensionless) and wall time (seconds).
    """
    n_modes = 2 * n_cells
    ops = build_ops(n_modes, 2)
    h = driven_hamiltonian(ops, CANDIDATES[sector], 0.5)
    rho0 = qt.ket2dm(initial_state(n_cells, 2, "edge_site"))
    baths = [DrudeLorentzBath(ops[j] + ops[j].dag(), lam=LAM,
                              gamma=GAMMA_CUTOFF, T=T_BATH, Nk=nk,
                              tag=f"site{j}")
             for j in range(n_modes)]
    t0 = time.time()
    solver = HEOMSolver(h, baths, max_depth=depth,
                        options={"nsteps": 20000, "progress_bar": False})
    result = solver.run(rho0, [0.0, PERIOD])
    state = result.states[-1]
    occ = np.asarray([float(np.real(qt.expect(op.dag() * op, state)))
                      for op in ops])
    assert np.all(np.isfinite(occ)), f"non-finite HEOM state {sector} d{depth}"
    assert abs(float(state.tr()) - 1.0) < 1e-6, "HEOM trace drift"
    return {
        "negativity_edge": float(negativity(state, 0, 31)),
        "negativity_center": float(negativity(state, n_modes // 2, 31)),
        "edge_population": float(occ[0] + occ[-1]),
        "n_ados": int(np.prod(solver.ados.dims)),
        "wall_seconds": time.time() - t0,
    }


def main() -> None:
    """Run the depth ladder for both sectors plus the Nk sensitivity point."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-cells", type=int, default=2, choices=(2,))
    p.add_argument("--depths", type=int, nargs="+", default=list(DEPTHS))
    p.add_argument("--output", type=Path,
                   default=ROOT / "data" / "heom_convergence.json")
    args = p.parse_args()

    start = time.time()
    points = {}
    for sector in CANDIDATES:
        for depth in args.depths:
            key = f"{sector}_depth{depth}_Nk{NK}"
            print(f"HEOM {key} ...", flush=True)
            points[key] = run_point(args.n_cells, sector, depth, NK)
            print(f"  neg_edge={points[key]['negativity_edge']:.6f} "
                  f"({points[key]['wall_seconds']:.0f}s)", flush=True)
    sens = {}
    for sector in CANDIDATES:
        key = f"{sector}_depth1_Nk{NK_SENSITIVITY}"
        print(f"HEOM {key} ...", flush=True)
        sens[key] = run_point(args.n_cells, sector, 1, NK_SENSITIVITY)
        print(f"  neg_edge={sens[key]['negativity_edge']:.6f} "
              f"({sens[key]['wall_seconds']:.0f}s)", flush=True)

    gaps = {}
    gates = {}
    for sector in CANDIDATES:
        n1 = points[f"{sector}_depth1_Nk{NK}"]["negativity_edge"]
        n2 = points[f"{sector}_depth2_Nk{NK}"]["negativity_edge"]
        n3 = points[f"{sector}_depth3_Nk{NK}"]["negativity_edge"]
        gaps[sector] = {"depth2_minus_depth1": n2 - n1,
                        "depth3_minus_depth2": n3 - n2}
        gates[f"{sector}_shrinking"] = bool(abs(n3 - n2) <= abs(n2 - n1))
        s1 = sens[f"{sector}_depth1_Nk{NK_SENSITIVITY}"]["negativity_edge"]
        gaps[sector]["Nk3_minus_Nk2_at_depth1"] = s1 - n1

    out = {
        "status": ("HEOM_CONVERGENCE_COMPLETE"
                   if all(gates.values()) else "HEOM_CONVERGENCE_FAILED"),
        "claim_scope": "hierarchy-depth convergence of the HEOM solver on "
                       "the production driven two-cell Kerr chain; solver "
                       "machinery validation only",
        "not_supported": ["topological protection theorem",
                          "production Lindblad re-validation",
                          "experimental validation", "metrological gain"],
        "design": {
            "n_cells": args.n_cells, "n_max": 2, "kerr": 0.5,
            "initial_state_mode": "edge_site",
            "bath": "Drude-Lorentz on every site, Q = a_j + a_j^dagger",
            "lam": LAM, "gamma_cutoff": GAMMA_CUTOFF, "T_bath": T_BATH,
            "Nk": NK, "Nk_sensitivity": NK_SENSITIVITY,
            "depths": list(args.depths), "t_stop": float(PERIOD),
            "disorder_level": 0.0,
        },
        "points": points,
        "matsubara_sensitivity": sens,
        "depth_gaps": gaps,
        "gates": gates,
        "software": {"qutip": qt.__version__, "numpy": np.__version__},
        "elapsed_seconds": time.time() - start,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(f"gates: {gates}")
    print(f"saved {args.output}")
    if not all(gates.values()):
        raise SystemExit("FAIL: " + ", ".join(k for k, v in gates.items()
                                              if not v))


if __name__ == "__main__":
    main()
