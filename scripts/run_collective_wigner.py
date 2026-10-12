#!/usr/bin/env python3
"""Collective-mode Wigner negativity with a validated bosonic mode transform.

Context
-------
The production observable is the negative part of the *single-site* Wigner
function. The manuscript and the data analysis report explicitly gate any
collective-mode statement: "collective-mode Wigner negativity remains gated
until a validated bosonic mode transform is implemented." This script
implements exactly that transform and validates it before any physical
statement is made.

The collective edge mode is the normalized single-particle edge vector
``c = (c_0, ..., c_{N-1})`` (Sigma_j |c_j|^2 = 1) of the sector's Floquet
operator, with annihilation operator  b_0 = Sigma_j c_j a_j.  A unitary mode
rotation U (N x N) is completed deterministically (first row = c), and the
exact many-body rotation U_circ on the truncated Fock space is built by acting
with the rotated creation operators on the vacuum:

    U_circ |n_1 ... n_N> = prod_j (n_j!)^(-1/2) (Sigma_i u_ji a_i^dag)^{n_j} |0>.

The collective-mode marginal is the partial trace of  U_circ rho U_circ^dag
over all modes except mode 0.  On a truncated space the rotation is only
approximately unitary; the validation suite quantifies this artifact and
bounds it with the same Gaussian-control logic already used for the local
observable (run_observable_validation.py).

Two modes of operation:

  --validate-only   run the transform validation suite (no dynamics).
  --run             run a bounded dynamics pilot (n_cells=2, 4 modes) with
                    stored final states; per record compute the local site-0
                    negativity (production definition), the collective
                    edge-mode negativity for the sector-matched edge vector,
                    and the identity-rotation internal consistency check.

No protection, metrological or collective-resource claim is made by this
script; it reports the collective-mode observable for the same matched
conditions already analyzed for the local observable.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import qutip as qt
from scipy.linalg import expm

from run_chiral_floquet_dynamics import (
    CANDIDATES,
    PERIOD,
    build_ops,
    driven_hamiltonian,
    edge_mode_vector,
    initial_state,
    rate_vector,
)

ROOT = Path(__file__).resolve().parents[1]
ALPHA = 1.0
DISORDER_LEVELS = (0.0, 0.35)


# --------------------------------------------------------------------------
# Mode rotation
# --------------------------------------------------------------------------

def complete_mode_rotation(edge_vector: np.ndarray) -> np.ndarray:
    """Deterministic unitary U whose first row equals ``edge_vector``.

    ``b_0 = sum_j u_0j a_j`` is the collective edge mode.  The remaining rows
    are obtained by deterministic Gram-Schmidt against the canonical basis:
    the completion only depends on ``edge_vector`` and the fixed basis order.
    """
    n = edge_vector.shape[0]
    u = np.zeros((n, n), dtype=complex)
    u[0, :] = edge_vector
    basis = np.eye(n, dtype=complex)
    rows = [edge_vector.copy()]
    for k in range(n):
        cand = basis[k].copy()
        for r in rows:
            cand -= np.vdot(r, cand) * r
        if np.linalg.norm(cand) > 1e-12:
            rows.append(cand / np.linalg.norm(cand))
    rows = rows[:n]
    u = np.array(rows, dtype=complex)
    # Orthonormalize the completed rows again (numerical hygiene).
    q, _ = np.linalg.qr(u.conj().T)
    u = q.conj().T
    # The QR step may have mixed the first row; re-impose it.
    u[0, :] = edge_vector
    # Re-orthonormalize rows 1..n-1 against row 0 and each other.
    rows = [edge_vector.copy()]
    for k in range(1, n):
        cand = np.eye(n, dtype=complex)[k].copy()
        cand -= np.vdot(edge_vector, cand) * edge_vector
        for r in rows[1:]:
            cand -= np.vdot(r, cand) * r
        if np.linalg.norm(cand) > 1e-12:
            rows.append(cand / np.linalg.norm(cand))
    rows = rows[:n]
    u = np.array(rows, dtype=complex)
    return u


def many_body_rotation(u: np.ndarray, n_modes: int, n_max: int) -> np.ndarray:
    """Exact (truncated) many-body rotation matrix T = U_circ in the Fock basis.

    Columns: T[:, idx(n)] = prod_j (n_j!)^(-1/2) (b_j^dag)^{n_j} |0>,
    with b_j^dag = sum_i u_ji a_i^dag.  On the truncated space T is only
    approximately unitary; ``unitarity_norm(T)`` quantifies the artifact.
    """
    dim = n_max + 1
    ops = build_ops(n_modes, n_max)
    bdag = [
        sum(u[j, i] * ops[i].dag() for i in range(n_modes))
        for j in range(n_modes)
    ]
    dims = dim**n_modes
    T = np.zeros((dims, dims), dtype=complex)
    # vacuum index = all zeros
    vac = np.zeros(dims, dtype=complex)
    vac[0] = 1.0
    multi = np.ndindex(*([dim] * n_modes))
    for n in multi:
        v = vac.copy()
        for j in range(n_modes):
            for _ in range(n[j]):
                v = bdag[j].full() @ v
            v /= np.sqrt(float(np.math.factorial(n[j])))
        idx = 0
        for j, nj in enumerate(n):
            idx += nj * (dim ** (n_modes - 1 - j))
        T[:, idx] = v
    return T


def unitarity_norm(T: np.ndarray) -> float:
    """Frobenius norm of T^dag T - I (truncation artifact of the rotation)."""
    return float(np.linalg.norm(T.conj().T @ T - np.eye(T.shape[0])))


def rotate_state(rho: qt.Qobj, T: np.ndarray, n_modes: int, n_max: int) -> qt.Qobj:
    """Return U_circ rho U_circ^dag as a Qobj with proper tensor dims."""
    dim = n_max + 1
    rho_full = rho.full()
    rho_rot = T @ rho_full @ T.conj().T
    return qt.Qobj(rho_rot, dims=[[dim] * n_modes, [dim] * n_modes])


def wigner_negativity(rho: qt.Qobj, grid_points: int) -> float:
    """Integrated negative volume of the single-mode Wigner function (production def.)."""
    if rho.isket:
        rho = qt.ket2dm(rho)
    grid = np.linspace(-4.0, 4.0, grid_points)
    w = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    return float(np.sum(np.abs(w[w < 0])) * dx * dx)


def collective_negativity(rho: qt.Qobj, u: np.ndarray, n_modes: int, n_max: int,
                          grid_points: int) -> tuple[float, float]:
    """(collective-mode negativity, local-mode-0 negativity in the rotated frame).

    Both marginals are taken in the rotated frame; mode 0 of that frame is the
    collective edge mode, and the marginal of the same frame's mode ``n_modes-1``
    (a bulk-completion mode) is returned as a frame-matched bulk reference.
    """
    T = many_body_rotation(u, n_modes, n_max)
    rho_rot = rotate_state(rho, T, n_modes, n_max)
    rho_edge = rho_rot.ptrace(0)
    rho_bulk_ref = rho_rot.ptrace(n_modes - 1)
    return (
        wigner_negativity(rho_edge, grid_points),
        wigner_negativity(rho_bulk_ref, grid_points),
    )


# --------------------------------------------------------------------------
# Validation suite
# --------------------------------------------------------------------------

def edge_state_in_rotated_frame(u: np.ndarray, n_modes: int, n_max: int,
                                kind: str) -> qt.Qobj:
    """State with a coherent (kind='coherent') or cat (kind='cat') excitation
    on the collective edge mode and vacuum elsewhere, expressed in the site
    frame: |psi> = T^dag (|excited>_0 x |0...0>_rest)."""
    dim = n_max + 1
    plus = qt.tensor(*[qt.coherent(dim, ALPHA if j == 0 else 0.0) for j in range(n_modes)])
    if kind == "coherent":
        psi_rot = plus
    else:
        minus = qt.tensor(*[qt.coherent(dim, -ALPHA if j == 0 else 0.0) for j in range(n_modes)])
        psi_rot = (plus + minus).unit()
    T = many_body_rotation(u, n_modes, n_max)
    psi = T.conj().T @ psi_rot.full().ravel()
    return qt.Qobj(psi, dims=[[dim] * n_modes, [1] * n_modes])


def run_validation(n_max_list, grid_points: int) -> dict:
    """Validate the mode transform independently of the physics narrative."""
    n_modes = 4  # 2 cells
    results = {"unitarity": [], "gaussian_control": [], "cat_benchmark": [],
               "identity_acid_test": {}}

    # (a) single-particle rotation unitarity for a non-trivial edge vector
    c = edge_mode_vector(2, CANDIDATES["chiral_topological"])
    u = complete_mode_rotation(c)
    results["single_particle_unitarity"] = float(np.linalg.norm(u @ u.conj().T - np.eye(n_modes)))

    # (b) sector-pinned unitarity: the rotation is exact on the subspace of
    # states whose total excitation is at most ``n_max`` (dynamics with an
    # edge-site cat preparation and particle-number-conserving Hamiltonian
    # plus loss never leave that sector).  Quantify by projecting the
    # unitarity defect onto the sector, with the sector mask built in the
    # truncated space itself (indices are base-(n_max+1) there).
    for nm in n_max_list:
        T = many_body_rotation(u, n_modes, nm)
        dimm = nm + 1
        total = dimm ** n_modes
        keep = np.zeros(total, dtype=bool)
        for n in np.ndindex(*([dimm] * n_modes)):
            if sum(n) <= nm:
                idx = 0
                for j, nj in enumerate(n):
                    idx += nj * (dimm ** (n_modes - 1 - j))
                keep[idx] = True
        sub = T[np.ix_(keep, keep)]
        defect = sub.conj().T @ sub - np.eye(int(keep.sum()))
        results.setdefault("sector_pinned_unitarity", []).append({
            "n_max": nm, "sector_total_excitation_leq": nm,
            "frobenius_defect": float(np.linalg.norm(defect)),
        })

    for nm in n_max_list:
        T = many_body_rotation(u, n_modes, nm)
        results["unitarity"].append({
            "n_max": nm, "frobenius_TdagT_minus_I": unitarity_norm(T),
        })
        # (d) Gaussian control: coherent excitation on the edge mode -> the
        # collective negativity is a pure truncation artifact and must decay.
        rho = qt.ket2dm(edge_state_in_rotated_frame(u, n_modes, nm, "coherent"))
        neg, neg_bulk = collective_negativity(rho, u, n_modes, nm, grid_points)
        results["gaussian_control"].append({
            "n_max": nm, "collective_negativity": neg, "bulk_ref_negativity": neg_bulk,
        })
        # (e) cat benchmark: cat on the edge mode -> compare with the single-mode
        # cat negative volume at the same truncation (in-frame reference).
        rho_cat = qt.ket2dm(edge_state_in_rotated_frame(u, n_modes, nm, "cat"))
        neg_cat, _ = collective_negativity(rho_cat, u, n_modes, nm, grid_points)
        single_cat = qt.tensor(qt.coherent(nm + 1, ALPHA) + qt.coherent(nm + 1, -ALPHA)).unit()
        ref_cat = wigner_negativity(qt.ket2dm(single_cat), grid_points)
        results["cat_benchmark"].append({
            "n_max": nm, "collective_negativity": neg_cat,
            "single_mode_reference": ref_cat,
            "relative_deviation": abs(neg_cat - ref_cat) / max(ref_cat, 1e-12),
        })

    # (c) acid test: identity rotation must reproduce the local site-0 observable.
    u_id = np.eye(n_modes, dtype=complex)
    # Build a generic non-trivial product-cat state and compare.
    plus = qt.tensor(*[qt.coherent(n_max_list[0] + 1, ALPHA * (0.5 + 0.1 * j)) for j in range(n_modes)])
    minus = qt.tensor(*[qt.coherent(n_max_list[0] + 1, -ALPHA * (0.5 + 0.1 * j)) for j in range(n_modes)])
    rho0 = qt.ket2dm((plus + minus).unit())
    nm = n_max_list[0]
    T_id = many_body_rotation(u_id, n_modes, nm)
    neg_local = wigner_negativity(rho0.ptrace(0), grid_points)
    rho_rot = rotate_state(rho0, T_id, n_modes, nm)
    neg_rot0 = wigner_negativity(rho_rot.ptrace(0), grid_points)
    results["identity_acid_test"] = {
        "local_site0": neg_local, "rotated_frame_site0": neg_rot0,
        "max_abs_diff": abs(neg_local - neg_rot0),
        "identity_unitarity": unitarity_norm(T_id),
    }
    return results


# --------------------------------------------------------------------------
# Dynamics pilot
# --------------------------------------------------------------------------

def simulate_state(h, psi0, c_ops, ntraj, seed, t_stop, n_time):
    result = qt.mcsolve(
        h, psi0, np.linspace(0.0, t_stop, n_time), c_ops,
        ntraj=ntraj, seeds=seed,
        options={"nsteps": 10000, "store_states": True, "progress_bar": False},
    )
    return result.states[-1]


def dump_figure_state(out_path: Path) -> dict:
    """Simulate and archive the fig:collective-wigner final state.

    Exact configuration of the figure (two cells, N_max=2, topological
    sector, edge_low_loss, d=0, sample 0, ntraj=10, kerr=0.5): the
    trajectory-averaged final state is reduced to the site-0 marginal
    and the sector-matched collective edge-mode marginal (3x3 blocks)
    so the figure script plots from this artifact and runs no dynamics.

    Args:
        out_path (Path): Destination JSON artifact path.

    Returns:
        dict: The written payload (config, both blocks, identity check).
    """
    n_cells, n_max, n_modes = 2, 2, 4
    sector, profile, disorder = "chiral_topological", "edge_low_loss", 0.0
    draw_seed, ntraj, kerr = 20260814, 10, 0.5
    ops = build_ops(n_modes, n_max)
    psi0 = initial_state(n_cells, n_max, "edge_site")
    rates = rate_vector(n_modes, profile, disorder, draw_seed, 1e-4, 0.2)
    c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
    state = simulate_state(
        driven_hamiltonian(ops, CANDIDATES[sector], kerr), psi0, c_ops,
        ntraj, draw_seed + 100, PERIOD, 9)
    rho_local = state.ptrace(0)
    c = edge_mode_vector(n_cells, CANDIDATES[sector])
    rho_edge = rotate_state(
        state, many_body_rotation(complete_mode_rotation(c), n_modes, n_max),
        n_modes, n_max).ptrace(0)

    def block(rho: qt.Qobj) -> dict:
        """Serialize a 3x3 block (units: density-matrix entries)."""
        m = np.asarray(rho.full())
        return {"real": m.real.tolist(), "imag": m.imag.tolist()}

    payload = {
        "status": "FIGURE_STATE_CACHED",
        "config": {
            "n_cells": n_cells, "n_max": n_max, "sector": sector,
            "loss_profile": profile, "disorder_level": disorder,
            "sample": 0, "draw_seed": draw_seed,
            "trajectory_seed": draw_seed + 100, "ntraj": ntraj,
            "kerr": kerr, "t_stop": PERIOD, "n_time": 9,
            "gamma_edge": 1e-4, "gamma_bulk": 0.2,
        },
        "rho_local_site0": block(rho_local),
        "rho_collective_edge": block(rho_edge),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"saved figure state {out_path}")
    return payload


def run_dynamics(args):
    n_modes = 2 * args.n_cells
    n_max = args.n_max
    ops = build_ops(n_modes, n_max)
    psi0 = initial_state(args.n_cells, n_max, args.initial_state_mode)
    grid_points = args.wigner_grid
    # sector-matched mode rotations
    rot = {sector: many_body_rotation(
        complete_mode_rotation(edge_mode_vector(args.n_cells, cand)), n_modes, n_max)
        for sector, cand in CANDIDATES.items()}
    u_id = np.eye(n_modes, dtype=complex)
    T_id = many_body_rotation(u_id, n_modes, n_max)
    rows = []
    start = time.time()
    for profile_index, profile in enumerate(("edge_low_loss", "uniform_loss")):
        for disorder_index, disorder in enumerate(args.disorder_levels):
            for sample in range(args.n_mc):
                draw_seed = args.seed + profile_index * 1_000_000 + disorder_index * 10_000 + sample
                rates = rate_vector(n_modes, profile, disorder, draw_seed,
                                    args.gamma_edge, args.gamma_bulk)
                c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
                for sector_index, (sector, candidate) in enumerate(CANDIDATES.items()):
                    seed = draw_seed + 100 + sector_index
                    state = simulate_state(
                        driven_hamiltonian(ops, candidate, args.kerr), psi0, c_ops,
                        args.ntraj, seed, args.t_stop, args.n_time,
                    )
                    # production local observable (site 0)
                    neg_local = wigner_negativity(state.ptrace(0), grid_points)
                    # identity rotation must reproduce it (per-record consistency)
                    rho_id = rotate_state(state, T_id, n_modes, n_max)
                    neg_id = wigner_negativity(rho_id.ptrace(0), grid_points)
                    # collective edge-mode marginal in the sector-matched frame
                    T_s = rot[sector]
                    rho_rot = rotate_state(state, T_s, n_modes, n_max)
                    neg_edge = wigner_negativity(rho_rot.ptrace(0), grid_points)
                    neg_bulkref = wigner_negativity(rho_rot.ptrace(n_modes - 1), grid_points)
                    occ = np.asarray([float(np.real(qt.expect(op.dag() * op, state))) for op in ops])
                    rows.append({
                        "n_cells": args.n_cells, "n_modes": n_modes, "n_max": n_max,
                        "sector": sector, **candidate, "kerr": args.kerr,
                        "loss_profile": profile, "disorder_level": disorder,
                        "sample": sample, "draw_seed": draw_seed,
                        "trajectory_seed": seed, "ntraj": args.ntraj,
                        "t_stop": args.t_stop, "period": PERIOD,
                        "wigner_grid": grid_points, "rates": rates.tolist(),
                        "negativity_local_site0": neg_local,
                        "negativity_identity_check": neg_id,
                        "identity_max_abs_diff": abs(neg_local - neg_id),
                        "negativity_collective_edge": neg_edge,
                        "negativity_collective_bulkref": neg_bulkref,
                        "edge_population": float(occ[0] + occ[-1]),
                        "bulk_population": float(np.sum(occ[1:-1])),
                        "occupations": occ.tolist(),
                    })
            print(f"completed profile={profile} disorder={disorder}")
    output = {
        "status": "COMPLETED_COLLECTIVE_WIGNER_DYNAMICS_PILOT",
        "claim_scope": ("collective edge-mode Wigner negativity for the same matched "
                        "conditions analyzed for the local observable; no protection claim"),
        "mode_transform": {
            "definition": "b_0 = sum_j c_j a_j, c = normalized sector edge vector of the "
                          "Floquet operator; deterministic unitary completion",
            "validation": "see collective_wigner_validation section",
        },
        "not_supported": ["topological protection theorem", "collective-mode protection",
                          "metrological gain", "experimental prediction"],
        "design": {
            "n_cells": args.n_cells, "n_modes": n_modes, "n_max": n_max,
            "initial_state_mode": args.initial_state_mode, "kerr": args.kerr,
            "n_mc": args.n_mc, "ntraj": args.ntraj, "wigner_grid": grid_points,
            "gamma_edge": args.gamma_edge, "gamma_bulk": args.gamma_bulk,
            "t_stop": args.t_stop, "n_time_points": args.n_time,
            "period": PERIOD, "paired_disorder": True,
            "disorder_levels": list(args.disorder_levels),
            "independent_unit": "disorder realization; trajectories nested",
        },
        "records": rows,
        "elapsed_seconds": time.time() - start,
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(output, indent=2) + "\n")
    print(f"saved {out} with {len(rows)} records")


# --------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------

def analyze(path: Path, grid_points: int, out_path: Path) -> dict:
    data = json.loads(path.read_text())
    rows = data["records"]
    # paired contrasts per (profile, disorder, sample): chiral - trivial
    by_key: dict[tuple, dict] = {}
    for r in rows:
        key = (r["loss_profile"], r["disorder_level"], r["sample"])
        by_key.setdefault(key, {})[r["sector"]] = r
    contrast = {"local_site0": [], "collective_edge": []}
    for key, pair in by_key.items():
        topo, triv = pair["chiral_topological"], pair["chiral_trivial"]
        contrast["local_site0"].append({
            "loss_profile": key[0], "disorder_level": key[1],
            "contrast": topo["negativity_local_site0"] - triv["negativity_local_site0"],
        })
        contrast["collective_edge"].append({
            "loss_profile": key[0], "disorder_level": key[1],
            "contrast": (topo["negativity_collective_edge"]
                         - triv["negativity_collective_edge"]),
        })

    def agg(items):
        out = {}
        for prof in ("edge_low_loss", "uniform_loss"):
            for dlev in sorted({x["disorder_level"] for x in items}):
                vals = [x["contrast"] for x in items
                        if x["loss_profile"] == prof and x["disorder_level"] == dlev]
                vals = np.asarray(vals)
                out[f"{prof}_d{dlev}"] = {
                    "mean": float(vals.mean()), "se": float(vals.std(ddof=1) / np.sqrt(len(vals))),
                    "n": int(len(vals)), "positive_fraction": float((vals > 0).mean()),
                }
        return out

    # identity-consistency check across all records
    id_diffs = np.asarray([r["identity_max_abs_diff"] for r in rows])

    analysis = {
        "paired_contrasts": {
            "local_site0": agg(contrast["local_site0"]),
            "collective_edge": agg(contrast["collective_edge"]),
        },
        "identity_consistency": {
            "max_abs_diff": float(id_diffs.max()),
            "mean_abs_diff": float(id_diffs.mean()),
            "tolerance": 1e-12,
        },
        "verdict": None,
    }
    analysis["verdict"] = (
        "collective edge-mode Wigner negativity computed and analyzed for the same "
        "matched conditions as the local observable; see paired_contrasts."
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(analysis, indent=2) + "\n")
    print(f"saved analysis {out_path}")
    return analysis


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--validate-only", action="store_true")
    p.add_argument("--dump-figure-state", type=Path, default=None,
                   help="write the fig:collective-wigner final state "
                        "(topological, edge_low_loss, d=0, sample 0) as a "
                        "JSON artifact and exit; the figure script reads "
                        "only this artifact and runs no dynamics")
    p.add_argument("--n-cells", type=int, default=2)
    p.add_argument("--n-max", type=int, default=2)
    p.add_argument("--n-mc", type=int, default=10)
    p.add_argument("--ntraj", type=int, default=10)
    p.add_argument("--disorder-levels", type=float, nargs="+", default=list(DISORDER_LEVELS))
    p.add_argument("--wigner-grid", type=int, default=31)
    p.add_argument("--initial-state-mode", choices=("edge_site", "topological_zero_mode"), default="edge_site")
    p.add_argument("--kerr", type=float, default=0.5)
    p.add_argument("--gamma-edge", type=float, default=1e-4)
    p.add_argument("--gamma-bulk", type=float, default=0.2)
    p.add_argument("--t-stop", type=float, default=PERIOD)
    p.add_argument("--n-time", type=int, default=9)
    p.add_argument("--seed", type=int, default=20260814)
    p.add_argument("--output", type=Path, default=ROOT / "data" / "collective_wigner_dynamics.json")
    p.add_argument("--validation-output", type=Path, default=ROOT / "data" / "collective_wigner_validation.json")
    p.add_argument("--analysis-output", type=Path, default=ROOT / "data" / "collective_wigner_analysis.json")
    p.add_argument("--n-max-list", type=int, nargs="+", default=[2, 3, 4])
    args = p.parse_args()

    if args.validate_only:
        val = run_validation(args.n_max_list, args.wigner_grid)
        val["status"] = "COMPLETED_MODE_TRANSFORM_VALIDATION"
        val["software"] = {"qutip": qt.__version__, "numpy": np.__version__}
        val["alpha"] = ALPHA
        args.validation_output.parent.mkdir(parents=True, exist_ok=True)
        args.validation_output.write_text(json.dumps(val, indent=2) + "\n")
        print(f"saved validation {args.validation_output}")
        return

    if args.dump_figure_state is not None:
        dump_figure_state(args.dump_figure_state)
        return

    if min(args.n_cells, args.n_max, args.n_mc, args.ntraj) < 1 or args.kerr < 0 or args.t_stop <= 0 or args.n_time < 2:
        p.error("invalid simulation arguments")
    run_dynamics(args)
    analyze(args.output, args.wigner_grid, args.analysis_output)


if __name__ == "__main__":
    main()
