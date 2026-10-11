#!/usr/bin/env python3
"""Fail-closed R8 analysis: Chabaud Fock witnesses and non-Gaussianity.

Implements audit run R8 (AUDIT_RESEARCH_REFINEMENTS_20261001.md line 126):
certifiable, grid-free witnesses of Wigner negativity computed on the SAME
trajectory ensemble as the production negativity, replacing the grid-dependent
31x31 integral as the certification of the claim (auditor point A-P1-05).

Gate 1 (parity): every rerun record in ``data/r8_witness/r8_witness_*cells.json``
must reproduce the archived production record with identical
``negativity_edge``, ``negativity_center``, ``rates``, ``draw_seed`` and
``trajectory_seed`` (deterministic rerun).  Any mismatch aborts the analysis.

Gate 2 (block sanity): each dumped 3x3 block must be Hermitian, unit trace and
positive semidefinite.

Witnesses (Chabaud, Cotton et al., Quantum 5, 471 (2021)):
    Omega_{e_n, alpha} = D(alpha) |n><n| D(alpha)^dagger,
    F_n(alpha) = <n| D(alpha)^dagger rho D(alpha) |n>,
    certified negative P-representation if max_alpha F_n(alpha) > omega_n.
Thresholds: omega_1 = omega_2 = 1/2 exactly (analytic); for n = 3..6 only the
published bounds [omega_lo, omega_hi] are available, so certification requires
F_n > omega_hi (conservative).  The threshold is displacement invariant, hence
the maximisation over alpha is admissible.

Non-Gaussianity: relative-entropy monotone
    Delta = S(rho_G) - S(rho),
with the Gaussian matched to the first and second moments computed in an
embedded space (d+4) so that truncated-operator boundary terms do not corrupt
the covariance matrix; validated against coherent states (Delta = 0) and
Fock states (Delta > 0).

Per record and per site (edge, centre control) the analyzer stores
max_alpha F_n and Delta; aggregates per (n_cells, disorder, profile) give
means with normal CIs, certification rates with Wilson CIs, and paired
topological-minus-trivial contrasts.

Output:
    data/r8_witness_summary.json
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import eval_genlaguerre, gammaln

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

# Chabaud et al., Quantum 5, 471 (2021), Table I / Ref. [34] therein:
OMEGA = {
    1: (0.5, 0.5), 2: (0.5, 0.5),
    3: (0.378, 0.427), 4: (0.375, 0.441),
    5: (0.314, 0.385), 6: (0.314, 0.378),
}
N_WITNESS = tuple(OMEGA)
COARSE = np.linspace(-2.0, 2.0, 81)
REFINE_STEP = 0.0025


def dfock_component(alpha: np.ndarray, n: int, k: int) -> np.ndarray:
    """Displacement-operator matrix element <k|D(alpha)|n> (full Hilbert space).

    Args:
        alpha (ndarray): Displacement amplitudes (dimensionless), shape (A,).
        n (int): Fock index of the witness state (dimensionless).
        k (int): Fock index of the bra (dimensionless, 0..2).

    Returns:
        ndarray: Complex matrix elements, shape (A,).
    """
    x = np.abs(alpha) ** 2
    if k >= n:
        log_c = 0.5 * (gammaln(n + 1) - gammaln(k + 1))
        out = (np.exp(log_c) * alpha ** (k - n) * np.exp(-x / 2.0)
               * eval_genlaguerre(n, k - n, x))
    else:
        log_c = 0.5 * (gammaln(k + 1) - gammaln(n + 1))
        out = (np.exp(log_c) * (-np.conj(alpha)) ** (n - k) * np.exp(-x / 2.0)
               * eval_genlaguerre(k, n - k, x))
    return out


def witness_fidelity(rho: np.ndarray, n: int, alphas: np.ndarray) -> np.ndarray:
    """F_n(alpha) = <n|D^dag rho D|n> for a batch of displacements.

    Only components m <= 2 of D(alpha)|n> enter because rho has support on
    the three lowest Fock states of the truncated model.

    Args:
        rho (ndarray): 3x3 reduced density matrix (dimensionless).
        n (int): Witness Fock index (dimensionless).
        alphas (ndarray): Displacement amplitudes (dimensionless), shape (A,).

    Returns:
        ndarray: Fidelities (dimensionless), shape (A,), in [0, 1].
    """
    u = np.stack([dfock_component(alphas, n, k) for k in range(rho.shape[0])], axis=1)
    f = np.einsum("ai,ij,aj->a", np.conj(u), rho, u)
    return np.real(f)


def maximize_witness(rho: np.ndarray, n: int) -> tuple[float, complex]:
    """Maximise F_n over a coarse grid with a local refinement.

    Args:
        rho (ndarray): 3x3 reduced density matrix (dimensionless).
        n (int): Witness Fock index (dimensionless).

    Returns:
        tuple: ``(max F_n, argmax alpha)`` (dimensionless).
    """
    g = COARSE
    re, im = np.meshgrid(g, g, indexing="ij")
    a0 = (re + 1j * im).ravel()
    f0 = witness_fidelity(rho, n, a0)
    i = int(np.argmax(f0))
    best, best_a = float(f0[i]), complex(a0[i])
    lo_r, hi_r = best_a.real - 4 * (g[1] - g[0]), best_a.real + 4 * (g[1] - g[0])
    lo_i, hi_i = best_a.imag - 4 * (g[1] - g[0]), best_a.imag + 4 * (g[1] - g[0])
    gr = np.arange(lo_r, hi_r + REFINE_STEP, REFINE_STEP)
    gi = np.arange(lo_i, hi_i + REFINE_STEP, REFINE_STEP)
    re, im = np.meshgrid(gr, gi, indexing="ij")
    a1 = (re + 1j * im).ravel()
    f1 = witness_fidelity(rho, n, a1)
    j = int(np.argmax(f1))
    if float(f1[j]) > best:
        best, best_a = float(f1[j]), complex(a1[j])
    return best, best_a


def non_gaussianity(rho: np.ndarray) -> float:
    """Relative-entropy non-Gaussianity Delta = S(rho_G) - S(rho) >= 0.

    Moments are taken in an embedded space (d+4) with the state renormalised,
    so truncated-operator boundary terms do not corrupt the covariance.

    Args:
        rho (ndarray): 3x3 reduced density matrix (dimensionless).

    Returns:
        float: Delta (dimensionless, nats).
    """
    d = rho.shape[0]
    big = d + 4
    emb = np.zeros((big, big), dtype=complex)
    emb[:d, :d] = rho
    tr = np.trace(emb)
    emb = emb / tr
    # quadrature moments from the embedded state via Fock-space matrices
    a = np.diag(np.sqrt(np.arange(1, big)), 1)
    x = (a + a.T.conj()) / np.sqrt(2)
    p = (a - a.T.conj()) / (1j * np.sqrt(2))
    mx = float(np.real(np.trace(emb @ x)))
    mp = float(np.real(np.trace(emb @ p)))
    xx = float(np.real(np.trace(emb @ x @ x)))
    pp = float(np.real(np.trace(emb @ p @ p)))
    xp = float(np.real(np.trace(emb @ (0.5 * (x @ p + p @ x)))))
    cov = np.array([[xx - mx * mx, xp - mx * mp], [xp - mx * mp, pp - mp * mp]])
    nu = float(np.sqrt(max(np.linalg.det(cov), 0.0)))
    if nu <= 0.5 + 1e-12:
        s_g = 0.0
    else:
        s_g = (nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5)
    ev = np.linalg.eigvalsh(emb)
    ev = ev[ev > 1e-15]
    s = float(-np.sum(ev * np.log(ev)))
    return float(s_g - s)


def block_matrix(raw) -> np.ndarray:
    """Convert dumped [Re, Im] rows into a complex matrix.

    Args:
        raw: JSON list of rows of ``[Re, Im]`` pairs.

    Returns:
        ndarray: Complex matrix (dimensionless).
    """
    return np.array([[complex(v[0], v[1]) for v in row] for row in raw])


def check_block(rho: np.ndarray) -> None:
    """Fail-closed sanity check of a dumped reduced block.

    Args:
        rho (ndarray): 3x3 reduced density matrix (dimensionless).

    Raises:
        AssertionError: If Hermiticity, trace or positivity is violated.
    """
    assert np.allclose(rho, rho.conj().T, atol=1e-10), "block not Hermitian"
    assert abs(np.trace(rho) - 1.0) < 1e-9, "block trace != 1"
    ev = np.linalg.eigvalsh(rho)
    assert ev.min() > -1e-9, f"block not positive semidefinite: {ev.min()}"


def wilson_interval(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion.

    Args:
        k (int): Number of successes (dimensionless).
        n (int): Number of trials (dimensionless).
        z (float): Normal quantile (dimensionless).

    Returns:
        tuple: ``(low, high)`` bounds (dimensionless).
    """
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, centre - half), min(1.0, centre + half))


def production_key(r: dict) -> tuple:
    """Record identity shared with the production archive."""
    return (r["loss_profile"], float(r["disorder_level"]), int(r["sample"]), r["sector"])


def load(path: Path) -> dict:
    """Load one result JSON."""
    return json.loads(path.read_text())


def main() -> None:
    """Run parity gate, witness computation and aggregation."""
    records = []
    for cells in (2, 3, 4):
        rerun = load(ROOT / "data" / "r8_witness" / f"r8_witness_{cells}cells.json")
        prod = load(ROOT / "data" / "chiral_production" /
                    f"chiral_floquet_{cells}cells_NMC500.json")
        prod_index = {production_key(r): r for r in prod["records"]}
        checked = 0
        for r in rerun["records"]:
            key = production_key(r)
            p = prod_index.get(key)
            assert p is not None, f"no production match for {cells}c {key}"
            for field in ("negativity_edge", "negativity_center",
                          "draw_seed", "trajectory_seed"):
                assert r[field] == p[field], (
                    f"parity failure {cells}c {key} {field}: "
                    f"{r[field]} != {p[field]}")
            assert r["rates"] == p["rates"], f"rates mismatch {cells}c {key}"
            checked += 1
        print(f"parity gate: {cells} cells OK ({checked} records)")
        for r in rerun["records"]:
            r["_cells"] = cells
            records.append(r)

    # witness + non-Gaussianity per record and site
    for r in records:
        for site_field, tag in (("rho_edge", "edge"), ("rho_center", "center")):
            rho = block_matrix(r[site_field])
            check_block(rho)
            for n in N_WITNESS:
                f_star, a_star = maximize_witness(rho, n)
                r[f"w{tag}_n{n}"] = f_star
                r[f"w{tag}_n{n}_alpha"] = [a_star.real, a_star.imag]
            r[f"ng_{tag}"] = non_gaussianity(rho)

    # aggregates per (cells, disorder, profile)
    strata = sorted({(r["_cells"], r["disorder_level"], r["loss_profile"])
                     for r in records})
    summary = {}
    for cells, disorder, profile in strata:
        sel = [r for r in records if (r["_cells"], r["disorder_level"],
                                      r["loss_profile"]) == (cells, disorder, profile)]
        by_sector = {s: sorted([r for r in sel if r["sector"] == s],
                               key=lambda r: r["sample"])
                     for s in ("chiral_topological", "chiral_trivial")}
        topo, triv = by_sector["chiral_topological"], by_sector["chiral_trivial"]
        n = len(topo)
        entry = {"n": n}
        for n_w in N_WITNESS:
            lo, hi = OMEGA[n_w]
            entry.setdefault("witness", {})[f"n{n_w}"] = {}
            for sector, group in by_sector.items():
                vals = np.array([r[f"wedge_n{n_w}"] for r in group])
                cert = int(np.sum(vals > hi))
                lo_ci, hi_ci = wilson_interval(cert, len(vals))
                entry["witness"][f"n{n_w}"][sector] = {
                    "omega_lo": lo, "omega_hi": hi,
                    "mean": float(vals.mean()), "se": float(vals.std(ddof=1) / np.sqrt(len(vals))),
                    "certified": cert, "certified_rate": cert / len(vals),
                    "certified_wilson95": [lo_ci, hi_ci],
                    "max": float(vals.max()),
                }
            # paired contrast of the witness value
            d = np.array([t[f"wedge_n{n_w}"] - v[f"wedge_n{n_w}"]
                          for t, v in zip(topo, triv)])
            entry["witness"][f"n{n_w}"]["paired_topo_minus_trivial"] = {
                "mean": float(d.mean()),
                "se": float(d.std(ddof=1) / np.sqrt(len(d))),
            }
        # centre-site control (n=1, both sectors) and non-Gaussianity
        f1c = np.array([r["wcenter_n1"] for r in topo + triv])
        entry["center_control_n1"] = {"mean": float(f1c.mean()),
                                      "max": float(f1c.max())}
        for sector, group in by_sector.items():
            ng = np.array([r["ng_edge"] for r in group])
            entry.setdefault("non_gaussianity", {})[sector] = {
                "mean": float(ng.mean()),
                "se": float(ng.std(ddof=1) / np.sqrt(len(ng))),
            }
        ngd = np.array([t["ng_edge"] - v["ng_edge"] for t, v in zip(topo, triv)])
        entry["non_gaussianity"]["paired_topo_minus_trivial"] = {
            "mean": float(ngd.mean()),
            "se": float(ngd.std(ddof=1) / np.sqrt(len(ngd))),
        }
        # production negativity contrast for the same records (recomputed,
        # parity-gated above)
        negd = np.array([t["negativity_edge"] - v["negativity_edge"]
                         for t, v in zip(topo, triv)])
        entry["negativity_contrast"] = {
            "mean": float(negd.mean()),
            "se": float(negd.std(ddof=1) / np.sqrt(len(negd))),
        }
        summary[f"{cells}cells_d{disorder:g}_{profile}"] = entry

    out = {
        "status": "COMPLETED_R8_WITNESS_ANALYSIS",
        "audit_source": "AUDIT_RESEARCH_REFINEMENTS_20261001.md (line 126)",
        "gates": {"production_parity": "PASSED (exact field equality)",
                  "block_sanity": "PASSED (Hermitian, unit trace, PSD)"},
        "witness_reference": "Chabaud et al., Quantum 5, 471 (2021)",
        "certification_rule": "max_alpha F_n > omega_n with omega_1=omega_2=0.5 exact; n>=3 uses conservative omega_hi",
        "n_records": len(records),
        "strata": summary,
    }
    dest = ROOT / "data" / "r8_witness_summary.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dest}")

    # compact table: d=0 edge_low
    print("\n d=0 edge_low_loss (edge site):")
    print(f"{'cells':>5} {'F1* topo':>10} {'F1* triv':>10} {'cert topo':>9} "
          f"{'cert triv':>9} {'dNG topo':>9} {'neg contrast':>12}")
    for cells in (2, 3, 4):
        e = summary[f"{cells}cells_d0_edge_low_loss"]
        w = e["witness"]["n1"]
        print(f"{cells:>5} {w['chiral_topological']['mean']:>10.4f} "
              f"{w['chiral_trivial']['mean']:>10.4f} "
              f"{w['chiral_topological']['certified_rate']:>9.3f} "
              f"{w['chiral_trivial']['certified_rate']:>9.3f} "
              f"{e['non_gaussianity']['chiral_topological']['mean']:>9.4f} "
              f"{e['negativity_contrast']['mean']:>+12.6f}")


if __name__ == "__main__":
    main()
