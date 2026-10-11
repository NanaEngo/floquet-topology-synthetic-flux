# Edge-mode pilot parameter audit — 2026-08-14

## Executive verdict

The edge-mode pilot parameters are **numerically coherent normalized parameters**, but
they are not an experimentally calibrated SI set and are not optimal parameters.
They were selected to test the topology–dissipation mechanism and numerical
conditioning, not by a constrained physical optimization.

```text
parameter consistency: PASS as a normalized model
SI calibration:       NOT COMPUTED
hardware feasibility: NOT ESTABLISHED
optimality:           NOT OPTIMIZED
```

## Parameters actually used

| Quantity | Pilot value | Interpretation | Physical status |
|---|---:|---|---|
| `J_A` topological | 0.1 | dimensionless intracell pulse coupling | uncalibrated |
| `J_B` topological | 1.2 | dimensionless intercell pulse coupling | uncalibrated |
| `J_A` trivial | 1.2 | matched control | uncalibrated |
| `J_B` trivial | 0.1 | matched control | uncalibrated |
| Period | `2*pi` | normalized drive period | no SI reference declared |
| `gamma_edge` | `1e-4` | normalized edge Lindblad rate | uncalibrated |
| `gamma_bulk` | `0.2` | normalized bulk Lindblad rate | mechanism contrast, not a measured damping |
| Kerr `chi` | `0` or `0.5` | passive control / nonlinear perturbation | uncalibrated |
| coherent amplitude | `ALPHA=1` | normalized cat-state amplitude | photon/resource mapping missing |
| Fock cutoff | `n_max=2` | numerical truncation | convergence required |
| disorder | `0, 0.5, 1.0, 1.5` | fractional perturbation of the rates | exploratory stress levels |
| Wigner grid | 31 points | numerical quadrature | convergence required |

## Why these values are not yet physical

The code does not declare a reference angular frequency \(\omega_{\rm ref}\), a
laser frequency, an input-output convention, a measured mechanical hopping, or
hardware bounds. Consequently, a normalized value cannot be converted uniquely to
Hz, watts, photon number, or device loss.

For illustration only, if the unit were \(\omega_{\rm ref}=2\pi\times1\) GHz,
then a normalized rate of `0.2` would correspond to approximately
\(2\pi\times200\) MHz, whereas `1e-4` would correspond to approximately
\(2\pi\times100\) kHz. The first value is not a plausible mechanical damping
interpretation for the cited optomechanical anchors; this demonstrates why the
pilot must not be presented as calibrated hardware dynamics.

Likewise, `chi=0.5` would be an extremely large nonlinear scale under a naive
1-GHz interpretation. It is retained only as a normalized nonlinear sensitivity
condition until a physical Kerr coefficient and photon-number convention are
specified.

The existing cross-project calibration audit identifies the unresolved quantities:

- reference frequency and angular/cyclic convention;
- phase-controlled mechanical hopping \(J_m\);
- input-drive amplitude and laser power mapping;
- device-specific loss rates;
- thermal occupation and detector noise;
- admissible hardware ranges.

See the parent calibration records in
`track_D_chaos/PARAMETER_AUDIT_20260814.md` and
`physical_parameter_manifest.json`.

## Why these values are not optimal

No optimization was run for this pilot. The selected points were chosen because

1. the independent chiral invariant is stable;
2. the zero and pi gaps are open;
3. the finite-chain edge diagnostic is present for the topological point;
4. the control point has the swapped hopping hierarchy;
5. the computation remains tractable.

This is a scientific design selection, not a Pareto optimum. An optimality claim
would require an explicit objective, admissible physical ranges, uncertainty model,
resource constraints and an out-of-sample validation of the selected point.

## Required calibration/optimization before physical claims

1. Declare \(\omega_{\rm ref}\) and all angular-frequency conventions.
2. Supply measured or source-backed SI ranges for \(J_A,J_B,\gamma_j,\chi\), drive
   amplitude and disorder.
3. Define the objective: edge-mode survival, Fisher information, robustness, or a
   multi-objective compromise.
4. Define constraints: stability, power, temperature, cutoff, fabrication and
   detector limits.
5. Run a constrained, uncertainty-aware Pareto study over those ranges.
6. Validate the selected region on independent parameter draws and larger Fock
   cutoffs.

Until these steps are complete, the manuscript wording must remain **normalized
mechanism pilot** or **model-level diagnostic**, never **physically optimal regime**.
