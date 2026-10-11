# Deterministic unitary covariance benchmark — 14 August 2026

## Purpose

The open-system edge-mode pilot showed that the conditional passive reference was not
validated. Before attributing that failure to loss, Kerr coupling or finite-size
physics, the kinematic observable itself was tested against the same coherent
single-particle drive used to define the Floquet edge vector.

The benchmark contains no Lindblad losses, Kerr term, disorder or Monte Carlo
sampling. It tests only whether

\[
u_e(t)=U_{\mathrm{coh}}(t)u_e(0),\qquad
b_e(t)=\sum_j u_{e,j}(t)^*a_j
\]

is covariant with the implemented micromotion.

## Provenance

- Script: `scripts/run_covariant_unitary_benchmark.py`.
- SLURM wrapper: `scripts/trackA_covariant_unitary_benchmark.sbatch`.
- Full run: SLURM job `15357`.
- Environment: `/home/nanaengo/miniforge3/envs/qom/bin/python`.
- Full output: `data/edge_mode_unitary_benchmark.json`.
- Local smoke output: `data/edge_mode_unitary_smoke.json`.
- Time grid: 257 equally spaced points over one period.
- Sizes: 2, 3, 4 and 6 unit cells (4, 6, 8 and 12 modes).

## Definitions and thresholds

For the edge vector and a complete orthonormal basis of vectors orthogonal to it, the
benchmark records

\[
p_e(t)=|u_e(t)^\dagger U_{\mathrm{coh}}(t)u_e(0)|^2,
\]

and

\[
\ell_{e\mathrm b}(t)=\max_{v\perp u_e(0)}
|u_e(t)^\dagger U_{\mathrm{coh}}(t)v|^2.
\]

The prespecified thresholds were:

| Quantity | Threshold |
|---|---:|
| Floquet eigenvector residual | \(10^{-10}\) |
| Unitarity residual | \(10^{-11}\) |
| Covariant population error \(1-p_e\) | \(10^{-11}\) |
| Edge–bulk overlap \(\ell_{e\mathrm b}\) | \(10^{-11}\) |

## Results

| Cells | Covariant population error (maximum) | Edge–bulk overlap (maximum) | Fixed-frame population minimum | Verdict |
|---:|---:|---:|---:|:---:|
| 2 | \(8.88\times10^{-16}\) | \(1.42\times10^{-31}\) | 0.895483 | PASS |
| 3 | \(1.33\times10^{-15}\) | \(7.72\times10^{-32}\) | 0.894713 | PASS |
| 4 | \(8.88\times10^{-16}\) | \(1.28\times10^{-31}\) | 0.894709 | PASS |
| 6 | \(8.88\times10^{-16}\) | \(8.60\times10^{-32}\) | 0.894709 | PASS |

The covariant observable therefore preserves the transported edge population to
machine precision in the coherent unitary model. The fixed laboratory-frame
observable exhibits micromotion oscillations and is not an equivalent replacement.

## Scientific interpretation

This result closes one previously unresolved **kinematic** gate:

```text
Floquet micromotion/operator correspondence       PASS
unitary edge–bulk orthogonality                    PASS
open-system passive survival reference             NOT VALIDATED
Kerr dynamics                                      NOT TESTED HERE
collective-mode Wigner negativity                  NOT VALIDATED
topological protection                             NOT ESTABLISHED
```

The benchmark does not show that losses preserve the edge mode. It shows only that
any failure of the passive reference in the open system cannot be attributed to an
incorrect fixed-frame operator or to a failure of unitary micromotion covariance.
Residual open-system edge–bulk transfer, nonlinearity, truncation and dissipation
remain to be separated in a future controlled study.
