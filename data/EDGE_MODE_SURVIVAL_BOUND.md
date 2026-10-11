# Conditional edge-mode survival bound

**Status:** conditional analytical framework; pilot validation not passed; not a theorem for the full Kerr–Lindblad model.
**Parameter status:** normalized, uncalibrated and not optimized; physical conversion is explicitly out of scope for this bound.

## 1. Scope

The certified topology belongs to the coherent, linear, periodically driven SSH
backbone. The open nonlinear dynamics are a perturbation of that backbone. We do
not assign a topological invariant to the full nonlinear Liouvillian unless a
separate Liouvillian invariant is derived.

Let the one-period single-particle Floquet operator have an isolated edge mode
\(u_e\), normalized by \(u_e^\dagger u_e=1\), and define

\[
b_e=\sum_j u_{e,j}^*a_j,
\qquad
P_e=b_e^\dagger|0\rangle\langle0|b_e
\]

at the single-particle level. The effective loss rate associated with site losses
\(\gamma_j\) is

\[
\gamma_e=\sum_j\gamma_j|u_{e,j}|^2.
\]

This quantity is an observable-dependent exposure measure, not itself a topological
invariant.

## 2. Exact passive benchmark

In the isolated-mode, passive limit (no Kerr, no edge–bulk coupling during the
observation window, and a fixed Floquet phase convention), the mode amplitude obeys

\[
\dot\beta_e(t)=\left[-i\varepsilon_e-\frac{\gamma_e}{2}\right]\beta_e(t),
\]

and therefore

\[
 n_e(t)=n_e(0)\exp(-\gamma_e t).
\]

The production code records this expression as `passive_prediction_edge_population`.
The completed passive pilot (job `15343`) tested the residual, but the cutoff panel (job `15344`) produced non-monotone maximum residuals of 0.03883, 0.02821 and 0.03557 for `N_max=2,3,4`. These values are not currently assignable to a declared numerical-error budget. The passive reference therefore remains a validation target, not an established bound for the driven finite chain.

The passive pilot must test the residual

\[
R_e(t)=n_e^{\mathrm{sim}}(t)-n_e(0)e^{-\gamma_e t}
\]

against integration, trajectory, Fock-cutoff and Wigner-grid errors. A small residual
is required before the Kerr extension is interpreted.

## 3. Finite-dimensional Wigner observable

For a local or validated collective-mode reduced state in a Fock cutoff \(d\),
define

\[
\mathcal N_W^{(d)}(\rho)=\frac12\left(\int|W_\rho^{(d)}(\alpha)|\,d^2\alpha-1\right).
\]

For two states in the same finite-dimensional representation,

\[
|\mathcal N_W^{(d)}(\rho)-\mathcal N_W^{(d)}(\sigma)|
\leq \frac12\|W_\rho^{(d)}-W_\sigma^{(d)}\|_1.
\]

A trace-norm bound can be used only after the cutoff-dependent map constant and the
Wigner-grid quadrature error are controlled. The present pilot therefore records
local Wigner negativity but does not silently call it a collective-mode Wigner
negativity.

## 4. Perturbative extension target

For a full model with Kerr strength \(\chi\), disorder \(V(t)\), and residual
edge–bulk coupling \(K_{eb}(t)\), the desired result is a finite-time inequality of
the schematic form

\[
\|\rho(t)-\rho_e(t)\|_1
\leq \epsilon_{\mathrm{trunc}}(d,t)
+\epsilon_{\mathrm{grid}}(h,t)
+\int_0^t \!\left[
 C_1\|K_{eb}(s)\|
+C_2\|V_{\mathrm{sb}}(s)\|
+C_3\chi\,n_{\max}
\right]ds,
\]

where \(V_{\mathrm{sb}}\) denotes symmetry-breaking disorder. The constants and
norms must be derived for the declared finite-dimensional representation; they are
not assumed universal. Combining this with the Wigner continuity inequality gives
a conditional bound on \(\mathcal N_W^{(d)}\), not an asymptotic protection law.

## 5. Required validation sequence

The deterministic unitary micromotion/edge--bulk covariance benchmark has passed for
2, 3, 4 and 6 cells. It validates the transported observable but does not include
loss and therefore does not validate the survival reference.

1. `kerr=0`, zero disorder: verify the passive exponential benchmark.
2. `kerr=0`, symmetry-preserving disorder: quantify the increase in \(\gamma_e\) and
   edge–bulk leakage.
3. `kerr=0`, symmetry-breaking disorder: test loss of the edge-mode diagnostic.
4. `kerr>0`: compare with the passive bound only within a declared perturbative
   regime; do not treat deviations as proof of protection.
5. Repeat for \(N_{\max}=2,3,4\), multiple Wigner grids and several time windows.
6. Compare edge-site and collective-mode observables under the same matched
   preparation and resource budget.

## 6. Claim boundary

Until these tests pass, the allowed statement is:

> The certified chiral Floquet model supports an edge-mode population diagnostic and
> a conditional passive-limit survival reference. The finite-time Kerr–Lindblad
> dynamics remain a numerical mechanism test.

The following statements are prohibited:

- “topology protects Wigner negativity”;
- “the bound proves asymptotic protection”;
- “the full Liouvillian is topological”;
- “the result implies metrological enhancement.”
