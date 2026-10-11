# Edge-mode pilot results — 14 August 2026

## Provenance

- Mechanistic pilot array: SLURM job `15343`, tasks 0–3.
- Cutoff-convergence array: SLURM job `15344`, tasks 0–5.
- Environment: `/home/nanaengo/miniforge3/envs/qom`.
- Executable: `scripts/run_chiral_edge_mode_pilot.py`.
- Auditors: `scripts/analyze_edge_mode_pilot.py` and `scripts/analyze_edge_mode_convergence.py`.
- Protocol: `data/EDGE_MODE_MECHANISM_PROTOCOL.md`.
- Bound framework: `data/EDGE_MODE_SURVIVAL_BOUND.md`.

Both arrays completed with empty stderr logs. The four pilot files contain 160
records each and 80 matched sector pairs. The six convergence files contain six
records and three matched sector pairs each. All JSON files passed structural,
finite-value and sector-pairing audits.

## Design

The pilot used two cells/four modes, ten disorder realizations per condition and
five nested quantum trajectories. It crossed:

- common edge-site versus topological edge-mode preparation;
- `chi=0` versus `chi=0.5`;
- `edge_low_loss` versus `uniform_loss`;
- disorder levels 0, 0.5, 1.0 and 1.5.

The convergence panel used the common edge-site preparation, three matched units,
`edge_low_loss`, and `N_max=2,3,4` in both the passive and Kerr regimes.

## Audited findings

For the common edge-site preparation at zero disorder and edge-low loss, the
collective edge-population contrast (topological minus trivial) was:

| Regime | Mean contrast | 95% interval | Interpretation |
|---|---:|---:|---|
| Passive, `chi=0` | 0.02698 | [0.02466, 0.02930] | positive descriptive contrast |
| Kerr, `chi=0.5` | 0.00962 | [−0.00166, 0.02091] | uncertain, not a robust positive effect |

For the topological edge-mode preparation at the same condition, the corresponding
contrasts were `7.98e-6` (passive; 95% interval [−4.53e-6, 2.05e-5]) and `4.03e-6`
(Kerr; 95% interval [1.79e-6, 6.27e-6]). These small values do not establish a
mechanism-level separation. Local site-0 negativity was exactly zero in the
reported zero-disorder pilot rows for the topological edge-mode preparation; it is
therefore not a substitute for collective-mode Wigner negativity.

The passive-reference residuals in the cutoff panel were:

| `N_max` | Maximum absolute residual |
|---:|---:|
| 2 | 0.03883 |
| 3 | 0.02821 |
| 4 | 0.03557 |

The residual is neither monotone in `N_max` nor demonstrably within a declared
numerical-error budget. The present data therefore do **not** validate the exact
isolated-mode exponential reference for the driven finite chain. Kerr runs were
correctly marked `comparison_only` and were not evaluated against the passive bound.

## Gate verdict

```text
Structural output audit             PASS
Paired design and finite values      PASS
Passive-bound validation             NOT PASSED
Nmax=2→3→4 convergence              NOT ESTABLISHED
Topology-by-loss mechanism           NOT ESTABLISHED
Collective-mode Wigner negativity    NOT COMPUTED/VALIDATED
Physical calibration/optimality      NOT COMPUTED
```

The pilot is a useful implementation and failure diagnostic. It does not support
claims of topological protection, a causal edge-loss mechanism, a Liouvillian
invariant, metrological gain, or experimental feasibility. The appropriate next
step is to repair the observable/model correspondence—specifically, derive the
Floquet-periodic edge-mode projector and its passive evolution or add an explicit
edge–bulk-decoupling control—before any larger production. Increasing Monte Carlo
sampling under the present unresolved bound is not justified.
