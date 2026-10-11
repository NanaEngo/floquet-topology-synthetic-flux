# Micromotion-covariant edge-mode pilot — 14 August 2026

## Purpose of the correction

The first edge-mode pilot used a fixed laboratory-frame operator and a coefficient
convention retained from the legacy implementation. That observable was not
sufficiently tied to the time-periodic Floquet eigenvector. The corrected pilot
uses

\[
b_e(t)=\sum_j u_{e,j}(t)^*a_j,
\qquad
u_e(t)=U_{\mathrm{coh}}(t)u_e(0),
\]

with the standard number-operator coefficients

\[
b_e^\dagger(t)b_e(t)=\sum_{ij}u_{e,i}(t)u_{e,j}(t)^*a_i^\dagger a_j.
\]

The effective passive exposure is integrated over the micromotion,

\[
\Gamma_e(t)=\int_0^t\sum_j\gamma_j|u_{e,j}(s)|^2\,ds,
\qquad
n_e^{\mathrm{ref}}(t)=n_e(0)e^{-\Gamma_e(t)}.
\]

This correction does not alter the earlier outputs. It creates a separately
provenanced pilot that tests the physically appropriate time-dependent observable.

## Provenance and completeness

- SLURM array: `15353`.
- Script: `scripts/run_chiral_edge_mode_pilot.py`.
- Wrapper: `scripts/trackA_edge_mode_covariant.sbatch`.
- Observable frame: `micromotion`.
- Operator convention: `standard`.
- Environment: `/home/nanaengo/miniforge3/envs/qom`.
- Four tasks completed with empty stderr logs.
- Four JSON files contain 160 records each and 80 matched sector pairs.
- Audit: `data/edge_mode_covariant_analysis.json`.

The local smoke test passed before submission with two records, finite time series,
standard operator convention and explicit micromotion metadata.

## Results at zero disorder and edge-low loss

| Preparation | Regime | Mean topological−trivial edge-population contrast | 95% interval |
|---|---|---:|---:|
| Common edge site | Passive, `chi=0` | 0.000388 | [0.000358, 0.000418] |
| Common edge site | Kerr, `chi=0.5` | 0.00317 | [−0.00869, 0.01502] |
| Topological edge mode | Passive, `chi=0` | −1.59×10⁻⁵ | [−1.82×10⁻⁵, −1.36×10⁻⁵] |
| Topological edge mode | Kerr, `chi=0.5` | −8.89×10⁻⁶ | [−1.16×10⁻⁵, −6.15×10⁻⁶] |

The common-site passive contrast is small relative to the earlier fixed-frame
contrast. The mode-resolved preparation gives a small negative contrast rather than
a positive survival advantage.

For the corrected passive reference, the maximum absolute residual was:

- common edge-site preparation: `0.03956`;
- topological edge-mode preparation: `0.38257`.

These residuals do not meet a declared numerical-error budget and are not evidence
that the finite driven chain realizes the isolated-mode exponential model.

## Scientific verdict

```text
Micromotion/operator correction        IMPLEMENTED and smoke-tested
Structural and finite-value audit      PASS
Passive survival bound                 NOT PASSED
Mode-preparation mechanism              NOT SUPPORTED by this pilot
Topology-by-loss interaction            NOT ESTABLISHED
Collective-mode Wigner negativity       NOT COMPUTED/VALIDATED
Physical calibration/optimality         NOT COMPUTED
```

The corrected observable removes a laboratory-frame ambiguity but does not rescue
the proposed protection mechanism. The remaining residual is consistent with
edge–bulk coupling and/or an incomplete Floquet-mode projection in the dissipative
finite-chain implementation.

The next valid action was a deterministic unitary micromotion test and an explicit
edge–bulk-decoupling benchmark. That benchmark subsequently passed for 2, 3, 4 and
6 cells, with maximum covariant population error $1.33\times10^{-15}$ and maximum
edge–bulk overlap $1.42\times10^{-31}$. This closes the kinematic observable gate,
but not the open-system passive gate. No topological-protection, metrological or
experimental claim is supported.

