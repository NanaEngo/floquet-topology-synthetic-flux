"""Exact master-equation cross-check of the four-period revival (audit A6).

Compares QuTiP ``mesolve`` (exact Lindblad integration) against ``mcsolve``
(Monte Carlo wave-function, ``ntraj = 200`` trajectories) on one fixed
two-cell production realization (edge-selective profile, ``d = 0.0``,
``kerr = 0.5``, seed 20260817) propagated to 4 drive periods
(``tlist`` of 33 points, ``31 x 31`` Wigner grid, site-0 negativity).

Writes the per-sector summary to
``data/revision_scans/multiperiod_mesolve_crosscheck.json`` (or
``sys.argv[1]``) and prints the records to stdout. Confirms the 4T
revival against the exact solver on this realization.
"""
import json
import sys
import time

import numpy as np
import qutip as qt

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))

from run_chiral_floquet_dynamics import (
    ALPHA,
    CANDIDATES,
    PERIOD,
    HALF_PERIOD,
    build_ops,
    driven_hamiltonian,
    initial_state,
    negativity,
    rate_vector,
)


def mesolve_crosscheck():
    n_cells = 2
    n_max = 2
    n_modes = 2 * n_cells
    ntraj = 200
    periods = 4
    t_stop = periods * PERIOD
    grid = 31
    kerr = 0.5
    seed = 20260817
    ops = build_ops(n_modes, n_max)
    psi0 = initial_state(n_cells, n_max, "edge_site")
    draw_seed = seed  # edge_low_loss, d=0.0, sample 0
    rates = rate_vector(n_modes, "edge_low_loss", 0.0, draw_seed, 1e-4, 0.2)
    c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]

    records = []
    for sector_index, (sector, candidate) in enumerate(CANDIDATES.items()):
        h = driven_hamiltonian(ops, candidate, kerr)
        tlist = np.linspace(0.0, t_stop, 33)

        exact = qt.mesolve(h, psi0, tlist, c_ops, options={"nsteps": 20000})
        neg_exact = negativity(exact.states[-1], 0, grid)

        traj_seed = draw_seed + 100 + sector_index
        mc = qt.mcsolve(
            h, psi0, tlist, c_ops, ntraj=ntraj, seeds=traj_seed,
            options={"nsteps": 20000, "store_states": True, "progress_bar": False},
        )
        neg_mc = negativity(mc.states[-1], 0, grid)

        occ_exact = np.asarray([float(np.real(qt.expect(op.dag() * op, exact.states[-1]))) for op in ops])
        records.append({
            "n_cells": n_cells, "n_modes": n_modes, "n_max": n_max,
            "sector": sector, "periods": periods, "ntraj": ntraj,
            "draw_seed": draw_seed, "trajectory_seed": traj_seed,
            "loss_profile": "edge_low_loss", "disorder_level": 0.0,
            "sample": 0, "wigner_grid": grid, "kerr": kerr,
            "negativity_mesolve": neg_exact,
            "negativity_mcsolve_ntraj%d" % ntraj: neg_mc,
            "diff": float(abs(neg_exact - neg_mc)),
            "edge_population_mesolve": float(occ_exact[0] + occ_exact[-1]),
        })

    out = {
        "status": "MULTIPERIOD_MESOLVE_CROSSCHECK",
        "claim_scope": "physicality check of the four-period revival; exact Lindblad vs Monte Carlo on one realization",
        "design": {
            "n_cells": 2, "n_max": 2, "edge_site": True, "kerr": 0.5,
            "gamma_edge": 1e-4, "gamma_bulk": 0.2, "periods": 4,
            "ntraj": ntraj, "wigner_grid": grid, "seed": seed,
        },
        "records": records,
        "software": {"numpy": np.__version__, "qutip": qt.__version__},
        "elapsed_seconds": round(time.time() - __import__("time").time(), 1),
    }
    return out


if __name__ == "__main__":
    import time as _t
    start = _t.time()
    out = mesolve_crosscheck()
    out["elapsed_seconds"] = round(_t.time() - start, 1)
    path = sys.argv[1] if len(sys.argv) > 1 else "data/revision_scans/multiperiod_mesolve_crosscheck.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out["records"], indent=1))
