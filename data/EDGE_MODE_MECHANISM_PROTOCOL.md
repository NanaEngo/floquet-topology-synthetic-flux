# Edge-mode mechanism protocol

**Version:** 2026-08-14, pre-pilot freeze
**Status:** pilot completed and structurally audited; scientific continuation gates not passed
**Parameter status:** normalized, uncalibrated and not optimized; see `data/EDGE_MODE_PARAMETER_AUDIT_20260814.md`

## Primary question

Does the certified chiral Floquet sector alter the finite-time survival of a
mode-resolved quantum resource through reduced exposure to bulk dissipation?

## Factors

| Factor | Levels | Role |
|---|---|---|
| Floquet sector | `chiral_topological` \((\nu_0,\nu_\pi)=(-1,0)\); `chiral_trivial` \((0,0)\) | primary contrast |
| Loss profile | `edge_low_loss`; `uniform_loss` | mechanistic interaction |
| Preparation | common edge-site cat; topological edge-mode cat | matched versus mode-resolved control |
| Kerr | 0; 0.5 | passive benchmark versus nonlinear extension |
| Disorder | 0, 0.5, 1.0, 1.5 | prespecified robustness levels |

The first pilot uses two cells/four modes, ten paired disorder realizations per
condition and five nested quantum trajectories. The same disorder draw, initial
state, resource amplitude, Fock cutoff, time grid and loss rates are used for the
two sectors within each paired unit. `J_A`, `J_B`, `gamma_edge`, `gamma_bulk`,
`chi` and `ALPHA` are dimensionless model coordinates; no SI or hardware claim is
attached to them.

## Primary and secondary outcomes

Primary outcome:

\[
\Delta n_e(t)=n_e^{\mathrm{topological}}(t)-n_e^{\mathrm{trivial}}(t),
\qquad
n_e=\langle b_e^\dagger b_e\rangle.
\]

Secondary outcomes:

- \(\gamma_e=\sum_j\gamma_j|u_{e,j}|^2\);
- edge and bulk loss currents;
- edge-mode population and bulk-population proxy;
- local site-0 and centre Wigner negativity;
- passive-limit residual relative to \(n_e(0)e^{-\gamma_e t}\);
- site-resolved occupations.

A validated collective-mode Wigner transform is a separate gate. It must not be
replaced by site-0 negativity without explicit qualification.

## Causal controls

The following controls are mandatory before a mechanism claim:

1. same physical edge-site preparation in both sectors;
2. same topological edge-mode preparation vector applied to both sectors;
3. edge-low-loss versus uniform loss;
4. passive \(\chi=0\) versus Kerr \(\chi=0.5\);
5. symmetry-preserving versus symmetry-breaking disorder in a later stage;
6. open versus periodic boundary conditions in a later stage.

## Analysis rule

Contrasts are computed within each matched disorder realization before aggregation.
Report the mean, standard deviation, 95% interval and positive fraction. Quantum
trajectories are nested within a disorder realization and are not treated as
independent disorder replicates. The pilot is descriptive; no p-value is promoted
as evidence of protection without a preregistered effect-size criterion.

## Deterministic kinematic control

Before interpreting the open-system residuals, the micromotion-dependent observable
was tested in the lossless single-particle model for 2, 3, 4 and 6 cells. The
covariant population error was below $1.33\times10^{-15}$ and the edge--bulk overlap
below $1.42\times10^{-31}$ across 257 times in one period. This control passes the
observable/frame correspondence but does not replace the open-system passive test.

## Continuation gates

The pilot may proceed to larger dissipative sizes only if all of the following are satisfied; the completed jobs `15343` and `15344`, followed by covariant job `15353`, did not satisfy the passive-bound and convergence requirements:

- passive-limit residual is below the declared numerical-error budget;
- edge-mode population and \(\gamma_e\) agree with the certified Floquet vector;
- the topology-by-loss interaction has a reproducible direction;
- the direction survives the common edge-site preparation;
- the result is not an artefact of `N_max=2`, Wigner grid or trajectory count;
- no contradiction appears between local and mode-resolved observables.

Failure of the gate means the protocol is retained as a negative mechanism result;
additional Monte Carlo samples under the same design are not authorized as a remedy.
