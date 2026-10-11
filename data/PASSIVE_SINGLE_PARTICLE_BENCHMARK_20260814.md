# Deterministic passive single-particle benchmark — 14 August 2026

## Purpose

The covariant unitary benchmark established that the edge observable is implemented
correctly in the coherent model. The next question was whether the passive survival
reference fails only after nonlinear/Fock-space dynamics or already in the simplest
open single-particle problem.

We therefore propagated one excitation under the same two-step Hamiltonian with
site-resolved Lindblad loss represented by the corresponding non-Hermitian
single-particle propagator. The calculation contains no Kerr term, Fock truncation,
quantum trajectories or Wigner transform.

## Reference and observables

For the coherent edge vector transported as \(u_e(t)=U_{\mathrm{coh}}(t)u_e(0)\),
the integrated-exposure reference is

\[
n_e^{\mathrm{ref}}(t)=\exp\left[-\int_0^t
\sum_j\gamma_j|u_{e,j}(s)|^2\,ds\right].
\]

The single-excitation open evolution gives a one-particle component
\(\rho_1(t)=K(t)|u_e(0)\rangle\langle u_e(0)|K(t)^\dagger\). We record the
covariant edge population, the residual from the reference, and the population in
the instantaneous orthogonal complement of the edge vector.

## Provenance

- Script: `scripts/run_passive_single_particle_benchmark.py`.
- Wrapper: `scripts/trackA_passive_single_particle.sbatch`.
- Full run: SLURM job `15358`.
- Environment: `/home/nanaengo/miniforge3/envs/qom/bin/python`.
- Output: `data/passive_single_particle_benchmark.json`.
- Sizes: 2, 3, 4 and 6 cells (4, 6, 8 and 12 modes).
- Time grid: 257 points over one drive period.
- Rates: \(\gamma_{\mathrm{edge}}=10^{-4}\), \(\gamma_{\mathrm{bulk}}=0.2\).

## Results

| Cells | Modes | Loss profile | Maximum reference residual | Maximum bulk population |
|---:|---:|---|---:|---:|
| 2 | 4 | Edge-selective | \(9.79\times10^{-4}\) | \(1.31\times10^{-3}\) |
| 3 | 6 | Edge-selective | \(1.94\times10^{-4}\) | \(2.39\times10^{-4}\) |
| 4 | 8 | Edge-selective | \(1.94\times10^{-4}\) | \(2.38\times10^{-4}\) |
| 6 | 12 | Edge-selective | \(1.94\times10^{-4}\) | \(2.38\times10^{-4}\) |
| 2–6 | 4–12 | Uniform | \(\leq1.94\times10^{-15}\) | \(\leq4.33\times10^{-16}\) |

The coherent propagator unitarity residual remained below \(1.9\times10^{-15}\) in
all cases.

## Scientific interpretation

Uniform loss is proportional to the identity in the single-particle sector. It
therefore commutes with the coherent drive and provides an exact control for the
integrated-exposure reference. The nonzero residual and complementary-space
population under edge-selective loss show that the loss profile does not remain
separable from the micromotion: the open dynamics generate a small edge–bulk
component even though the unitary covariant observable is exact.

This result explains why the passive reference can fail in the finite driven chain,
but it does not establish a topological protection effect. It is a model-level
mechanistic decomposition:

```text
unitary observable/frame correspondence       PASS
uniform-loss passive reference                 PASS
structured-loss isolated-mode reference       NOT EXACT
Kerr/Fock/trajectory effects                   NOT TESTED HERE
collective-mode Wigner negativity              NOT VALIDATED
topological protection                         NOT ESTABLISHED
```

The structured-loss result should be presented as a deterministic diagnostic, not as
an estimate of a universal loss law or a causal proof of resource preservation.
