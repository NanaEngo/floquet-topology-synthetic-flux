"""Validate the exact <=2-quanta subspace reduction used by run_8mode_nmc500_parallel.py.

Checks: operator-level invariance vs the full-space model (Hamiltonian block
diff, coupling nnz, initial-state leakage), end-to-end subspace timing, and
record-level equivalence against the archived full-space rows (3-cell
production JSON and 4-cell partial checkpoint, seed-identical; tolerance set
by integrator branch noise, far below the production standard error).

Usage: mamba run -n qom python scripts/validate_subspace_reduction.py
"""
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import qutip as qt
import scipy.sparse as sp

sys.path.insert(0, "/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608/scripts")
import run_8mode_nmc500_parallel as R

ROOT = Path("/home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_A_topo_2608")
PERIOD = float(R.PERIOD)
Q_MAX = 2  # even cat: site 0 in {0,2}; H number-conserving; loss number-reducing


def sub_basis(n_modes, n_max, q_max=Q_MAX):
    tuples, idx = [], []
    for occ in itertools.product(range(n_max + 1), repeat=n_modes):
        if sum(occ) <= q_max:
            tuples.append(occ)
            v = 0
            for n in occ:
                v = v * (n_max + 1) + n
            idx.append(v)
    return tuples, np.array(idx, dtype=np.int64)


def sub_destroy(n_modes, n_max, tuples, site):
    lookup = {t: i for i, t in enumerate(tuples)}
    rows, cols, data = [], [], []
    for j, t in enumerate(tuples):
        if t[site] > 0:
            t2 = list(t)
            amp = np.sqrt(t[site])
            t2[site] -= 1
            i = lookup[tuple(t2)]
            rows.append(i)
            cols.append(j)
            data.append(amp)
    mat = sp.csr_matrix((data, (rows, cols)), shape=(len(tuples), len(tuples)))
    return qt.Qobj(mat)


def build_sub_ops(n_modes, n_max, tuples):
    return [sub_destroy(n_modes, n_max, tuples, j) for j in range(n_modes)]


def sub_step_hamiltonian(sub_ops, intra, inter, kerr):
    h = R.second_quantize(R.finite_chain_matrix(len(sub_ops) // 2, intra, inter), sub_ops)
    h += kerr * (sub_ops[0].dag() * sub_ops[0].dag() * sub_ops[0] * sub_ops[0])
    return h


def embed_dm(sub_qobj, idx, n_modes, n_max):
    """Embed a sub-basis density matrix back into the full tensor-product space."""
    d = sub_qobj.data.as_scipy().tocoo()
    rows, cols, vals = idx[d.row], idx[d.col], d.data
    n_full = (n_max + 1) ** n_modes
    mat = sp.csr_matrix((vals, (rows, cols)), shape=(n_full, n_full))
    dims = [[n_max + 1] * n_modes, [n_max + 1] * n_modes]
    return qt.Qobj(mat, dims=dims)


def simulate_sub(profile, disorder, sample, sector_key, n_cells, n_max, n_mc,
                 ntraj, wigner_grid, gamma_edge, gamma_bulk, t_stop, n_time,
                 kerr, base_seed):
    """Subspace twin of run_8mode_nmc500_parallel.simulate_single."""
    n_modes = 2 * n_cells
    tuples, idx = sub_basis(n_modes, n_max)
    ops = build_sub_ops(n_modes, n_max, tuples)
    candidate = R.CANDIDATES[sector_key]

    # initial state: full-space construction (same code), projected + renormalized
    psi0_full = R.initial_state(n_cells, n_max, "edge_site").full().ravel()
    outside = np.linalg.norm(np.delete(psi0_full, idx))
    psi0 = qt.Qobj(psi0_full[idx])

    profile_list = list(R.DISORDER_LEVELS)
    profile_index = ("edge_low_loss", "uniform_loss").index(profile)
    disorder_index = profile_list.index(disorder)
    draw_seed = base_seed + profile_index * 1_000_000 + disorder_index * 10_000 + sample
    rates = R.rate_vector(n_modes, profile, disorder, draw_seed, gamma_edge, gamma_bulk)
    c_ops = [np.sqrt(x) * op for x, op in zip(rates, ops) if x > 0]
    sector_index = list(R.CANDIDATES.keys()).index(sector_key)
    seed = draw_seed + 100 + sector_index

    ha = sub_step_hamiltonian(ops, candidate["intra_A"], 0.0, kerr)
    hb = sub_step_hamiltonian(ops, 0.0, candidate["inter_B"], kerr)
    h = qt.QobjEvo([[ha, R.coeff_A], [hb, R.coeff_B]])
    tlist = np.linspace(0.0, t_stop, n_time)
    t0 = time.time()
    result = qt.mcsolve(h, psi0, tlist, c_ops, ntraj=ntraj, seeds=seed,
                        options={"nsteps": 10000, "store_states": True, "progress_bar": False})
    solve_t = time.time() - t0
    state = embed_dm(result.states[-1], idx, n_modes, n_max)
    occ = [float(np.real(qt.expect(op.dag() * op, result.states[-1]))) for op in ops]
    neg_edge = R.negativity(state, 0, wigner_grid)
    neg_center = R.negativity(state, n_modes // 2, wigner_grid)
    row = {"negativity_edge": neg_edge, "negativity_center": neg_center,
           "occupations": occ, "edge_population": occ[0] + occ[-1],
           "bulk_population": float(np.sum(occ[1:-1])),
           "solve_seconds": solve_t, "psi0_outside_norm": outside}
    return row, len(tuples), idx, ops, psi0


def check_invariance(n_cells, n_max, tuples, idx, ops):
    n_modes = 2 * n_cells
    cand = R.CANDIDATES["chiral_topological"]
    full_ops = R.build_ops(n_modes, n_max)
    h_full = R.step_hamiltonian(full_ops, cand["intra_A"], 0.0, 0.5)
    h_sub = sub_step_hamiltonian(ops, cand["intra_A"], 0.0, 0.5)
    csr = h_full.data.as_scipy().tocsr()
    outside = np.setdiff1d(np.arange(csr.shape[0]), idx)
    cross = csr[outside][:, idx]
    block = csr[idx][:, idx]
    sub_csr = h_sub.data.as_scipy().tocsr()
    block_diff = abs(block - sub_csr).max() if block.nnz or sub_csr.nnz else 0.0
    # loss operators must not map sub states outside
    a_cross = 0
    for i in range(n_modes):
        a_csr = full_ops[i].data.as_scipy().tocsr()
        a_cross += abs(a_csr[outside][:, idx]).max() if a_csr[outside][:, idx].nnz else 0.0
    return cross.nnz, float(block_diff), float(a_cross)


def main():
    # --- 1) invariance at 4 cells (the production size) ---
    n_cells, n_max = 4, 3
    n_modes = 2 * n_cells
    tuples, idx = sub_basis(n_modes, n_max)
    ops = build_sub_ops(n_modes, n_max, tuples)
    cross_nnz, block_diff, a_cross = check_invariance(n_cells, n_max, tuples, idx, ops)
    print(f"[inv] sub_dim={len(tuples)} (full {(n_max+1)**n_modes})  "
          f"H cross-nnz={cross_nnz}  block|full-sub|={block_diff:.3e}  "
          f"loss cross max={a_cross:.3e}", flush=True)
    assert cross_nnz == 0, "Hamiltonian leaves the <=2-quanta subspace"
    assert block_diff < 1e-12, "subspace block mismatch"
    assert a_cross < 1e-12, "loss operator leaves the subspace"

    # --- 2) sub protocol at 4 cells, full production settings, one record ---
    row4, dim4, *_ = simulate_sub("edge_low_loss", 0.0, 0, "chiral_topological",
                                  4, 3, 500, 10, 31, 1e-4, 0.2, PERIOD, 9, 0.5, 20260814)
    print(f"[4c] dim={dim4} solve={row4['solve_seconds']:.2f}s "
          f"neg_e={row4['negativity_edge']:.6f} neg_c={row4['negativity_center']:.6f} "
          f"outside={row4['psi0_outside_norm']:.3e}", flush=True)

    # --- 3) record-level equivalence vs archived full-space rows ---
    # Tolerance rationale: the operator-level reduction is EXACT (invariance
    # block diff 0), so record diffs here are integrator step-path/trajectory
    # branch noise at different matrix dimensions (same seeds). Production
    # aggregates have SE >= 1e-3; record-level noise must stay well below.
    def compare(records, label, limit_n=12):
        diffs = []
        checked = 0
        for rec in records:
            if not (rec["loss_profile"] == "edge_low_loss" and rec["disorder_level"] == 0.0):
                continue
            row, *_ = simulate_sub(rec["loss_profile"], rec["disorder_level"], rec["sample"],
                                   rec["sector"], n_cells_cmp, 3, 500, rec["ntraj"],
                                   rec["wigner_grid"], 1e-4, 0.2, rec["t_stop"], 9,
                                   rec["kerr"], 20260814)
            dn_e = abs(row["negativity_edge"] - rec["negativity_edge"])
            dn_c = abs(row["negativity_center"] - rec["negativity_center"])
            d_occ = np.max(np.abs(np.array(row["occupations"]) - np.array(rec["occupations"])))
            diffs.append((rec.get("idx", checked), dn_e, dn_c, d_occ))
            checked += 1
            if checked >= limit_n:
                break
        arr = np.array([[d[1], d[2], d[3]] for d in diffs])
        print(f"[{label}] {checked} records: per-record max neg_e={arr[:,0].max():.2e} "
              f"neg_c={arr[:,1].max():.2e} occ={arr[:,2].max():.2e} "
              f"(median neg_e={np.median(arr[:,0]):.2e})", flush=True)
        assert arr[:, 0].max() < 1e-4, f"{label} neg_edge mismatch {arr[:,0].max()}"
        assert arr[:, 1].max() < 1e-4, f"{label} neg_center mismatch {arr[:,1].max()}"
        assert arr[:, 2].max() < 1e-3, f"{label} occupations mismatch {arr[:,2].max()}"

    n_cells_cmp = 3
    arch3 = json.loads((ROOT / "data/chiral_production_nmax3/chiral_floquet_3cells_NMC500_nmax3.json").read_text())["records"]
    compare(arch3, "3c-arch")

    n_cells_cmp = 4
    partial = []
    pf = ROOT / "data/chiral_production_nmax3/chiral_floquet_4cells_NMC500_nmax3.json.partial.jsonl"
    fin = ROOT / "data/chiral_production_nmax3/chiral_floquet_4cells_NMC500_nmax3.json"
    if pf.exists():
        for line in pf.read_text().splitlines():
            line = line.strip()
            if line:
                rec = json.loads(line)
                r = dict(rec["row"])
                r["idx"] = rec["idx"]
                partial.append(r)
    elif fin.exists():
        # final file is in canonical task order; the first 148 rows are the
        # full-space-provenance records (written before the subspace patch)
        for i, rec in enumerate(json.loads(fin.read_text())["records"][:148]):
            r = dict(rec)
            r["idx"] = i
            partial.append(r)
    if len(partial) >= 6:
        compare(partial, "4c-fullspace")
    else:
        print(f"[4c-fullspace] only {len(partial)} rows, skipping comparison", flush=True)
    print("ALL CHECKS PASSED", flush=True)


if __name__ == "__main__":
    main()
