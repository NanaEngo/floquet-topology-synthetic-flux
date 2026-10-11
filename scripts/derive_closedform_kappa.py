#!/usr/bin/env python3
"""Derive and archive the closed-form edge--bulk rate constants (SI Eq. S19).

Computes, at zero disorder under the production edge-low-loss profile, for
two, three, and four cells:

  Lambda_1          strict single-particle lower bound (SI text after Eq. S21):
                    (1/T) int_0^T dt sum_{j != 0} gamma_j |<j|U(t)|0>|^2
                    with U(t) the single-particle two-step evolution.
  Lambda_2^coh      rate-weighted two-photon dwell WITH the production Kerr
                    term (k_Kerr = 0.5 on the edge site), bosonic i<=j
                    convention (equivalently the ordered-pair sum over
                    distinguishable propagation):
                    (1/T) int_0^T dt sum_{(i,j)!=(0,0)} (gamma_i+gamma_j)
                        |<i,j|U_2(t)|0,0>|^2
  Lambda_2^coh,lin  same integral with the Kerr term omitted (linear chain).
  gamma_eb, kappa   measured residual from data/participation_bound_verification.json
                    and the per-size prefactor kappa = gamma_eb / Lambda_2^coh.

Trivial-sector control (SI paragraph "Sector dependence of the prefactor"):
  Lambda_2^coh(triv), gamma_eff(triv), gamma_meas(triv) from the archived
  N_MC=500 production records at d=0 edge-low loss, and the implied
  kappa_eff(triv) = (gamma_meas - gamma_eff) / Lambda_2^coh(triv).

Self-checks reproduce SI tab:closedform-kappa; the script aborts on any
mismatch. Output: data/closedform_kappa_derivation.json
"""
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from run_chiral_floquet_dynamics import (
    CANDIDATES,
    HALF_PERIOD,
    PERIOD,
    edge_mode_vector,
    finite_chain_matrix,
    rate_vector,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "closedform_kappa_derivation.json"

GAMMA_EDGE, GAMMA_BULK = 1e-4, 0.2
KERR = 0.5
PROFILE = "edge_low_loss"
N_CELLS_LIST = (2, 3, 4)
N_GRID = 4001          # odd -> a grid point sits exactly at t = HALF_PERIOD
TOL_LAMBDA2 = 5e-5     # reproduce SI tab:closedform-kappa to its 4th decimal
TOL_LAMBDA1 = 2e-4     # converged Lambda_1 vs previously quoted 4th decimal

# SI tab:closedform-kappa (previously quoted; Lambda_2 entries must match
# exactly, Lambda_1 must match within TOL_LAMBDA1).
SI_TABLE = {
    "lambda1": {2: 0.0129, 3: 0.0524, 4: 0.0524},
    "lambda2_coh": {2: 0.0120, 3: 0.0950, 4: 0.0950},
    "lambda2_coh_lin": {2: 0.0256, 3: 0.1046, 4: 0.1046},
    "kappa": {2: 11.8, 3: 4.2, 4: 4.2},
}


def rates(n_cells: int) -> np.ndarray:
    """Production edge-low-loss rates at zero disorder (units: 1/period)."""
    return rate_vector(2 * n_cells, PROFILE, 0.0, 0, GAMMA_EDGE, GAMMA_BULK)


def lambda1(n_cells: int) -> float:
    """Single-particle lower bound Lambda_1 (1/period)."""
    cand = CANDIDATES["chiral_topological"]
    ha = finite_chain_matrix(n_cells, cand["intra_A"], 0.0)
    hb = finite_chain_matrix(n_cells, 0.0, cand["inter_B"])
    gamma = rates(n_cells)
    n_modes = 2 * n_cells
    ts = np.linspace(0.0, PERIOD, N_GRID)
    dt = ts[1] - ts[0]
    u_a = expm(-1j * ha * dt)
    u_b = expm(-1j * hb * dt)
    u_ah = expm(-1j * ha * dt / 2.0)
    u_bh = expm(-1j * hb * dt / 2.0)
    psi = np.zeros(n_modes, dtype=complex)
    psi[0] = 1.0
    weight = gamma[1:]          # sum over j != 0
    total = 0.0
    for k in range(len(ts) - 1):
        step_a = ts[k] < HALF_PERIOD
        psi_mid = (u_ah if step_a else u_bh) @ psi
        total += float(np.sum(weight * np.abs(psi_mid[1:]) ** 2)) * dt
        psi = (u_a if step_a else u_b) @ psi
    return total / PERIOD


def two_photon_dwell(n_cells: int, candidate: dict, kerr: float) -> float:
    """Rate-weighted two-photon dwell Lambda_2 (1/period), ordered-pair sum."""
    n_modes = 2 * n_cells
    ha = finite_chain_matrix(n_cells, candidate["intra_A"], 0.0)
    hb = finite_chain_matrix(n_cells, 0.0, candidate["inter_B"])
    gamma = rates(n_cells)
    eye = np.eye(n_modes)
    h_a = np.kron(ha, eye) + np.kron(eye, ha)
    h_b = np.kron(hb, eye) + np.kron(eye, hb)
    # Kerr: kerr * n_0(n_0 - 1) acts only on the doubly occupied edge site.
    h_a[0, 0] += 2.0 * kerr
    h_b[0, 0] += 2.0 * kerr
    ts = np.linspace(0.0, PERIOD, N_GRID)
    dt = ts[1] - ts[0]
    u_a = expm(-1j * h_a * dt)
    u_b = expm(-1j * h_b * dt)
    u_ah = expm(-1j * h_a * dt / 2.0)
    u_bh = expm(-1j * h_b * dt / 2.0)
    psi = np.zeros(n_modes * n_modes, dtype=complex)
    psi[0] = 1.0                       # |0,0>: two quanta on the edge site
    idx, weight = [], []
    for i in range(n_modes):
        for j in range(n_modes):
            if i == 0 and j == 0:
                continue
            idx.append(i * n_modes + j)
            weight.append(gamma[i] + gamma[j])
    idx = np.array(idx)
    weight = np.array(weight)
    total = 0.0
    for k in range(len(ts) - 1):
        step_a = ts[k] < HALF_PERIOD
        psi_mid = (u_ah if step_a else u_bh) @ psi
        total += float(np.sum(weight * np.abs(psi_mid[idx]) ** 2)) * dt
        psi = (u_a if step_a else u_b) @ psi
    return total / PERIOD


def load_participation() -> dict:
    """Per-size topological-sector production rates (Gamma_eff, Gamma_meas)."""
    bound = json.loads(
        (ROOT / "data" / "participation_bound_verification.json").read_text()
    )["bound_verification"]
    return {row["n_cells"]: row for row in bound}


def trivial_decay(n_cells: int) -> dict:
    """Trivial-sector Gamma_eff and Gamma_meas at d=0, edge-low loss.

    Gamma_meas = -ln(mean Neg(T) / Neg(0)) / T from the archived N_MC=500
    production records; Neg(0) is the shared initial-state reference stored in
    the participation-bound artifact.
    """
    path = ROOT / "data" / "chiral_production" / f"chiral_floquet_{n_cells}cells_NMC500.json"
    records = json.loads(path.read_text())["records"]
    vals = [
        r["negativity_edge"] for r in records
        if r["sector"] == "chiral_trivial"
        and r["loss_profile"] == PROFILE
        and r["disorder_level"] == 0.0
    ]
    assert len(vals) == 500, f"expected 500 trivial d=0 records at {n_cells} cells, got {len(vals)}"
    neg0 = load_participation()[n_cells]["negativity_0"]
    gamma_meas = -np.log(float(np.mean(vals)) / neg0) / PERIOD
    gamma = rates(n_cells)
    u = edge_mode_vector(n_cells, CANDIDATES["chiral_trivial"])
    gamma_eff = float(np.sum(gamma * np.abs(u) ** 2))
    return {"gamma_eff": gamma_eff, "gamma_meas": float(gamma_meas), "n_records": len(vals)}


def main() -> None:
    participation = load_participation()
    sizes = {}
    failures = []

    for n in N_CELLS_LIST:
        l1 = lambda1(n)
        l2 = two_photon_dwell(n, CANDIDATES["chiral_topological"], KERR)
        l2_lin = two_photon_dwell(n, CANDIDATES["chiral_topological"], 0.0)
        row = participation[n]
        gamma_eff = float(row["gamma_eff"])
        gamma_meas = float(row["gamma_meas"])
        gamma_eb = gamma_meas - gamma_eff
        kappa = gamma_eb / l2

        triv = trivial_decay(n)
        l2_triv = two_photon_dwell(n, CANDIDATES["chiral_trivial"], KERR)
        kappa_eff_triv = (triv["gamma_meas"] - triv["gamma_eff"]) / l2_triv

        sizes[n] = {
            "lambda1": l1,
            "lambda2_coh": l2,
            "lambda2_coh_lin": l2_lin,
            "gamma_eff": gamma_eff,
            "gamma_meas": gamma_meas,
            "gamma_eb": gamma_eb,
            "kappa": kappa,
            "trivial": {
                "lambda2_coh": l2_triv,
                "gamma_eff": triv["gamma_eff"],
                "gamma_meas": triv["gamma_meas"],
                "gamma_eb": triv["gamma_meas"] - triv["gamma_eff"],
                "kappa_eff": kappa_eff_triv,
            },
        }

        # ---- self-checks against SI tab:closedform-kappa ----
        if abs(l2 - SI_TABLE["lambda2_coh"][n]) > TOL_LAMBDA2:
            failures.append(f"{n}c Lambda_2_coh {l2:.6f} != {SI_TABLE['lambda2_coh'][n]}")
        if abs(l2_lin - SI_TABLE["lambda2_coh_lin"][n]) > TOL_LAMBDA2:
            failures.append(f"{n}c Lambda_2_lin {l2_lin:.6f} != {SI_TABLE['lambda2_coh_lin'][n]}")
        if abs(l1 - SI_TABLE["lambda1"][n]) > TOL_LAMBDA1:
            failures.append(f"{n}c Lambda_1 {l1:.6f} != {SI_TABLE['lambda1'][n]} (tol {TOL_LAMBDA1})")
        if round(kappa, 1) != SI_TABLE["kappa"][n]:
            failures.append(f"{n}c kappa {kappa:.3f} != {SI_TABLE['kappa'][n]}")
        if not l1 < gamma_eb:
            failures.append(f"{n}c Lambda_1 {l1:.6f} !< gamma_eb {gamma_eb:.6f}")

    if abs(sizes[3]["lambda2_coh"] - sizes[4]["lambda2_coh"]) > 1e-6:
        failures.append("three/four-cell Lambda_2 degeneracy broken")
    if abs(sizes[3]["lambda1"] - sizes[4]["lambda1"]) > 1e-9:
        failures.append("three/four-cell Lambda_1 degeneracy broken")

    # SI quotes these measured decay rates (sector-dependence paragraph and
    # tab:closedform-kappa): they must reproduce from the production records.
    quoted_meas = {"topo": {2: 0.337, 3: 0.591, 4: 0.591},
                   "triv": {2: 0.532, 3: 0.533, 4: 0.533}}
    for n in N_CELLS_LIST:
        if round(sizes[n]["gamma_meas"], 3) != quoted_meas["topo"][n]:
            failures.append(f"{n}c topological gamma_meas {sizes[n]['gamma_meas']:.6f} != {quoted_meas['topo'][n]}")
        if round(sizes[n]["trivial"]["gamma_meas"], 3) != quoted_meas["triv"][n]:
            failures.append(f"{n}c trivial gamma_meas {sizes[n]['trivial']['gamma_meas']:.6f} != {quoted_meas['triv'][n]}")

    if failures:
        print("SELF-CHECK FAILURES:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)

    out = {
        "status": "COMPLETED_CLOSEDFORM_KAPPA_DERIVATION",
        "definitions": {
            "lambda1": "(1/T) int_0^T dt sum_{j!=0} gamma_j |<j|U(t)|0>|^2 (single-particle two-step evolution)",
            "lambda2_coh": "(1/T) int_0^T dt sum_{(i,j)!=(0,0)} (gamma_i+gamma_j)|<i,j|U_2(t)|0,0>|^2, bosonic i<=j convention, Kerr k=0.5 on the edge site",
            "lambda2_coh_lin": "same as lambda2_coh with the Kerr term omitted",
            "kappa": "gamma_eb / lambda2_coh, per-size prefactor (NOT size-universal)",
            "gamma_eb": "Gamma_meas - Gamma_eff from data/participation_bound_verification.json",
        },
        "model": {
            "n_cells": list(N_CELLS_LIST),
            "sector_topological": CANDIDATES["chiral_topological"],
            "sector_trivial": CANDIDATES["chiral_trivial"],
            "kerr": KERR,
            "profile": PROFILE,
            "gamma_edge": GAMMA_EDGE,
            "gamma_bulk": GAMMA_BULK,
            "disorder": 0.0,
            "period": PERIOD,
            "n_grid": N_GRID,
        },
        "sizes": {str(n): sizes[n] for n in N_CELLS_LIST},
        "previously_quoted": {
            "source": "SI tab:closedform-kappa (literals, no derivation archive)",
            "table": {k: {str(n): v[n] for n in N_CELLS_LIST} for k, v in SI_TABLE.items()},
            "lambda1_note": (
                "Converged recomputation gives 0.0128/0.0523/0.0523; the previously "
                "quoted 0.0129/0.0524/0.0524 differs by ~1e-4 (4th-decimal grid drift "
                "of the unarchived original computation). SI prose (0.013, 0.052) is "
                "unaffected. The SI table now carries the derived values."
            ),
            "sector_paragraph_note": (
                "The SI 'Sector dependence of the prefactor' paragraph quoted the "
                "withdrawn universal-kappa era values kappa~2.9, Lambda_2=0.05/0.24, "
                "kappa_eff=1.88/1.78/1.72 (equal to gamma_eb/Lambda_2 with the "
                "unreproducible old dwell integrals, e.g. 0.452560/0.2407=1.886 at "
                "two cells). Those constants are replaced by the derived values in "
                "this file: kappa_eff(triv)=2.86/2.69/2.48, "
                "Lambda_2(triv)=0.158/0.161/0.161; measured rates quoted in the "
                "paragraph (0.337/0.532/0.591/0.533) already matched production and "
                "are unchanged."
            ),
        },
        "self_checks": "all passed (Lambda_2 and Lambda_2_lin to 5e-5, Lambda_1 to 2e-4, kappa to 0.1, Lambda_1 < gamma_eb, 3/4-cell degeneracy, measured topological/trivial decay rates reproduce the SI-quoted 0.337/0.591 and 0.532/0.533)",
        "software": {"python": platform.python_version(), "numpy": np.__version__},
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    for n in N_CELLS_LIST:
        s = sizes[n]
        print(
            f"{n}c: L1={s['lambda1']:.6f} L2={s['lambda2_coh']:.6f} "
            f"L2lin={s['lambda2_coh_lin']:.6f} gamma_eb={s['gamma_eb']:.6f} "
            f"kappa={s['kappa']:.3f} | triv L2={s['trivial']['lambda2_coh']:.6f} "
            f"kappa_eff={s['trivial']['kappa_eff']:.3f}"
        )
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
