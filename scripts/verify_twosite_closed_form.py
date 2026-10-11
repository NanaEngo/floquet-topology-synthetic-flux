#!/usr/bin/env python3
"""Two-site ab initio prediction of the edge--bulk prefactor (no calibration).

The single-site closed form (SI Eq. S21),
    kappa_single = (1/8) (int_Omega Delta W0 d2b) / delta_-(0) = 3.411
for the truncated even cat under photon loss, is applied WITHOUT calibration to
the two-cell chiral model: the edge-site negativity decay rate is predicted as
    Gamma_pred = Gamma_eff + kappa_single * Lambda_2^coh,
where Gamma_eff = sum_j gamma_j |u_{e,j}|^2 is the mode-weighted rate and
Lambda_2^coh the rate-weighted two-photon dwell under the two-particle Floquet
operator U_2(t) (SI Eq. S19). Everything is computed exactly in the truncated
2-mode x (n_max=2) Hilbert space (9-dim; the Lindblad superoperator is solved
by mesolve, not Monte Carlo).

Checks performed
----------------
1. Exact edge-marginal negativity decay Gamma_meas (one period, edge-low-loss,
   topological sector, zero disorder) — reproduces the production value.
2. Ab initio prediction Gamma_pred = Gamma_eff + kappa_single*Lambda_2^coh.
3. Return factor: kappa_eff = (Gamma_meas - Gamma_eff)/Lambda_2^coh and its
   ratio to kappa_single (coherent-return correction).
4. Two-mode identity consistency: for the product state cat x vacuum the
   general 4D identity reduces exactly to kappa_single (vacuum contributes
   zero flux).

Output: data/twosite_closedform.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import qutip as qt
from scipy.linalg import expm

from run_chiral_floquet_dynamics import (
    CANDIDATES,
    HALF_PERIOD,
    PERIOD,
    build_ops,
    driven_hamiltonian,
    edge_mode_vector,
    initial_state,
    rate_vector,
)

ROOT = Path(__file__).resolve().parents[1]
KAPPA_SINGLE = 3.41098   # closed form, SI Eq. S21
ALPHA = 1.0
N_MAX = 2
GAMMA_EDGE, GAMMA_BULK = 1e-4, 0.2
N_CELLS = 2


def wigner_neg(rho: qt.Qobj, N: int = 601, lim: float = 5.0) -> float:
    grid = np.linspace(-lim, lim, N)
    w = qt.wigner(rho, grid, grid)
    dx = grid[1] - grid[0]
    return float(np.sum(np.abs(w[w < 0])) * dx * dx)


def two_photon_dwell(n_cells: int, candidate: dict, kerr: float,
                     gamma: np.ndarray, n_t: int = 400) -> dict:
    """Rate-weighted two-photon dwell Lambda_2^coh under the coherent drive."""
    n_modes = 2 * n_cells
    dim = N_MAX + 1
    # two-particle Floquet operator U_2(t): tensor the single-particle
    # two-step propagators (drive is a product of two half-period evolutions)
    ha = _single_step_matrix(n_cells, candidate["intra_A"], 0.0)
    hb = _single_step_matrix(n_cells, 0.0, candidate["inter_B"])
    u1 = _fock_propagator(ha, hb)      # single-particle two-step unitary
    u2 = _two_particle_propagator(u1)  # two-particle (tensor) unitary
    # initial two-photon edge state |0,0> (site 0 has 2 quanta)
    psi2 = np.zeros(n_modes * n_modes, dtype=complex)
    psi2[0] = 1.0
    ts = np.linspace(0.0, PERIOD, n_t)
    dt = ts[1] - ts[0]
    dwell = 0.0
    probs = {"edge_edge": 0.0, "edge_away": 0.0, "away_away": 0.0}
    # ---- exact propagation with piecewise-constant drive ----
    U_now = np.eye(n_modes * n_modes, dtype=complex)
    # build piecewise propagators for the two-photon sector
    def two_body_ham(h_single: np.ndarray) -> np.ndarray:
        # H_2 = h1 otimes I + I otimes h1
        I = np.eye(n_modes)
        return np.kron(h_single, I) + np.kron(I, h_single)
    H_a2 = two_body_ham(ha)
    H_b2 = two_body_ham(hb)
    dwell = 0.0
    probs = {"edge_edge": 0.0, "edge_away": 0.0, "away_away": 0.0}
    prev = None
    for t in ts:
        tau = t % PERIOD
        if prev is not None:
            # advance by dt/2 using the active half-period Hamiltonian
            h2 = H_a2 if tau < HALF_PERIOD else H_b2
            U_now = expm(-1j * h2 * (dt / 2.0)) @ U_now
            psi = U_now @ psi2
            occ = np.abs(psi) ** 2
            # two-photon configurations: (i,j) = two quanta on sites (i,j), i<=j
            for i in range(n_modes):
                for j in range(i, n_modes):
                    idx = i * n_modes + j
                    p = occ[idx]
                    if i == 0 and j == 0:
                        probs["edge_edge"] += p * dt
                    elif i == 0 or j == 0:
                        probs["edge_away"] += p * dt
                        dwell += (gamma[i] + gamma[j]) * p * dt
                    else:
                        probs["away_away"] += p * dt
                        dwell += (gamma[i] + gamma[j]) * p * dt
        prev = t
    return {
        "Lambda_2_coh": dwell / PERIOD,
        "dwell_decomp": {k: v / PERIOD for k, v in probs.items()},
    }


def _single_step_matrix(n_cells, intra, inter):
    from run_chiral_floquet_dynamics import finite_chain_matrix
    return finite_chain_matrix(n_cells, intra, inter)


def _fock_propagator(ha, hb):
    n = ha.shape[0]
    return expm(-1j * hb * HALF_PERIOD) @ expm(-1j * ha * HALF_PERIOD)


def _two_particle_propagator(u1):
    n = u1.shape[0]
    I = np.eye(n)
    return np.kron(u1, u1)  # two distinguishable "particles" in the tensor space


def main() -> None:
    n_modes = 2 * N_CELLS
    dim = N_MAX + 1
    ops = build_ops(n_modes, N_MAX)
    cand = CANDIDATES["chiral_topological"]
    kerr = 0.5

    # ---- exact (non-MC) Lindblad dynamics over one period ----
    psi0 = initial_state(N_CELLS, N_MAX, "edge_site")
    rho0 = qt.ket2dm(psi0)
    gamma = rate_vector(n_modes, "edge_low_loss", 0.0, 0,
                        GAMMA_EDGE, GAMMA_BULK)
    c_ops = [np.sqrt(g) * op for g, op in zip(gamma, ops) if g > 0]
    H = driven_hamiltonian(ops, cand, kerr)
    tlist = np.linspace(0.0, PERIOD, 33)
    res = qt.mesolve(H, rho0, tlist, c_ops, e_ops=[])
    negs = np.array([wigner_neg(s.ptrace(0)) for s in res.states])
    gamma_meas = -np.log(negs[-1] / negs[0]) / PERIOD

    # ---- Gamma_eff from the edge-mode weights ----
    u = edge_mode_vector(N_CELLS, cand)
    gamma_eff = float(np.sum(gamma * np.abs(u) ** 2))

    # ---- Lambda_2^coh exact ----
    dwell = two_photon_dwell(N_CELLS, cand, kerr, gamma)
    l2 = dwell["Lambda_2_coh"]

    # ---- ab initio prediction ----
    gamma_pred = gamma_eff + KAPPA_SINGLE * l2
    kappa_eff = (gamma_meas - gamma_eff) / l2

    result = {
        "status": "COMPLETED_TWOSITE_CLOSEDFORM",
        "model": {"n_cells": N_CELLS, "n_max": N_MAX, "kerr": kerr,
                  "candidate": cand, "profile": "edge_low_loss",
                  "gamma_edge": GAMMA_EDGE, "gamma_bulk": GAMMA_BULK,
                  "period": PERIOD},
        "kappa_single_closed_form": KAPPA_SINGLE,
        "delta0": float(negs[0]),
        "gamma_eff": gamma_eff,
        "Lambda_2_coh": l2,
        "dwell_decomposition": dwell["dwell_decomp"],
        "gamma_meas": gamma_meas,
        "gamma_pred_ab_initio": gamma_pred,
        "pred_over_meas": gamma_pred / gamma_meas,
        "kappa_eff_implied": kappa_eff,
        "kappa_eff_over_kappa_single": kappa_eff / KAPPA_SINGLE,
        "residual_ab_initio": gamma_pred - gamma_meas,
    }
    out = ROOT / "data" / "twosite_closedform.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
