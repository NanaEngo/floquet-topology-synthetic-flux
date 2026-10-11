# Independent chiral Floquet model manifest

**Status:** coherent bulk invariants and finite-chain edge diagnostics validated as model diagnostics; production dynamics audited with no size-robust local Wigner-protection signature.
**Parameter status:** normalized, uncalibrated and not optimized; see `data/EDGE_MODE_PARAMETER_AUDIT_20260814.md`.

## Model

The two-step drive has period `T=2*pi`:

- Step A, duration `pi`: intracell hopping `J_A`;
- Step B, duration `pi`: intercell hopping `J_B`;
- no onsite mass term;
- chiral operator `Gamma=sigma_z`.

The bulk Hamiltonian is

\[
H(k;J_A,J_B)=
\begin{pmatrix}
0 & J_A+J_B e^{-ik}\\
J_A+J_B e^{ik} & 0
\end{pmatrix},
\]

so `{Gamma,H}=0` exactly.

## Invariants

Symmetric time frames A and B are constructed from half-step propagators. The
winding of the off-diagonal block in each frame is `w_A` and `w_B`; the reported
integer invariants are

\[
\nu_0=(w_A+w_B)/2,\qquad
\nu_\pi=(w_A-w_B)/2.
\]

The implementation checks resolution 40, 80, 120 and 200 k-points, both Floquet
gaps, and the chiral-unitary residual in both symmetric frames.

## Selected points

| Sector | `J_A` | `J_B` | `nu_0` | `nu_pi` | minimum zero gap | minimum pi gap |
|---|---:|---:|---:|---:|---:|---:|
| Topological candidate | 0.1 | 1.2 | -1 | 0 | 2.1991 | 0.3142 |
| Trivial candidate | 1.2 | 0.1 | 0 | 0 | 2.1991 | 0.3142 |

The maximum chiral-frame residual is below `1e-15` for both selected points.
Finite-chain diagnostics at 4, 8 and 12 cells show two near-zero edge modes for
the topological candidate and no corresponding persistent near-zero pair for the
control candidate at larger size. This is consistent with, but is not by itself a
proof of, bulk--edge correspondence for the open nonlinear dynamics.

## Dissipative dynamics gate

A QuTiP pilot (job `15336`) used the same two sectors, paired disorder, edge-low-loss
and uniform-loss controls, 10 disorder realizations and 10 nested trajectories at 4
and 6 modes. The edge-site Wigner-negativity advantage was positive in several
4-mode edge-loss conditions but was not uniform at 6 modes and disappeared under
uniform loss. A Kerr-free sensitivity check also did not produce a robust positive
separation.

A bounded production array (job `15338`) then evaluated 2, 3 and 4 cells (4, 6
and 8 modes). The first two sizes used 500 paired disorder realizations; the 8-mode
case used 100 because of the QuTiP cost. All three outputs passed structural audits.
The edge-loss contrast was positive for 2 cells, negative for 3 cells, and negative
for 4 cells. Uniform-loss contrasts remained near zero or weakly negative.

Therefore:

```text
bulk topology: PASS
finite-chain edge diagnostic: PASS as a model diagnostic
production structural audit: PASS
size-robust local Wigner-protection signature: NOT SUPPORTED
large topological-protection claim: BLOCKED
```

These results are not connected to the archived Wigner arrays. No manuscript claim
has been upgraded beyond finite-time descriptive dynamics.

## Edge-mode mechanism follow-up

Because the local site-0 observable is not size-robust, the next pilot measures the
collective mode

```text
b_edge = sum_j conjugate(u_j) a_j
```

and its effective loss exposure

```text
gamma_edge_eff = sum_j gamma_j |u_j|^2.
```

The frozen design crosses common edge-site versus topological edge-mode cat
preparation with Kerr-free versus Kerr dynamics, paired edge-low-loss and
uniform-loss controls, and four prespecified disorder levels. The implementation is
`scripts/run_chiral_edge_mode_pilot.py`; the SLURM wrapper is
`scripts/trackA_edge_mode_pilot.sbatch`. The protocol and conditional passive bound
are `data/EDGE_MODE_MECHANISM_PROTOCOL.md` and
`data/EDGE_MODE_SURVIVAL_BOUND.md`.

The pilot array `15343` and cutoff-convergence array `15344` completed with empty
error logs and complete structural audits. The four pilot outputs each contain 160
records; the six convergence outputs each contain six records. The passive-reference
residuals were 0.03883, 0.02821 and 0.03557 at `N_max=2,3,4`, respectively. Because
these residuals were not monotone and no declared numerical-error budget was met,
the passive-bound and cutoff-convergence gates remain **NOT PASSED**. The pilot
results are summarized in `data/EDGE_MODE_PILOT_RESULTS_20260814.md` and audited in
`data/edge_mode_pilot_analysis.json` and
`data/edge_mode_convergence_analysis.json`.

Collective-mode Wigner negativity remains a separate validation gate; local site
negativity is not used as its substitute. No larger production or protection claim
is authorized until the edge-mode projector/evolution correspondence is repaired
and the passive benchmark is revalidated.

A corrected micromotion-covariant pilot was then executed as job `15353`. It uses the
standard `b_edge^dagger b_edge` convention and transports the edge vector through the
coherent Floquet micromotion. All four tasks completed with empty stderr logs. The
common-site passive contrast at zero disorder and edge-low loss was 0.000388, whereas
the mode-resolved passive contrast was −1.59e-5. Maximum passive residuals were
0.03956 and 0.38257, so the bound remains unvalidated. The detailed result and claim
boundary are in `data/EDGE_MODE_COVARIANT_PILOT_RESULTS_20260814.md`.

A deterministic unitary covariance benchmark was then completed in job `15357` for
2, 3, 4 and 6 cells. The covariant edge-population error remained below
$1.33\times10^{-15}$ and the transported edge--bulk overlap below
$1.42\times10^{-31}$. This validates the kinematic observable/frame correspondence
only; the open-system passive reference remains unvalidated. Full provenance is in
`data/EDGE_MODE_UNITARY_BENCHMARK_20260814.md`.

A deterministic one-excitation passive benchmark (job `15358`) isolates the loss
profile. Uniform loss reproduces the integrated-exposure reference to
$1.94\times10^{-15}$, whereas edge-selective loss leaves a residual of
$9.79\times10^{-4}$ at 4 modes and $1.94\times10^{-4}$ at 6--12 modes, with a
corresponding complementary-space population. This is a finite-model
loss--micromotion diagnostic, not a protection result. The full record is
`data/PASSIVE_SINGLE_PARTICLE_BENCHMARK_20260814.md`.
