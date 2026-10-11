# Floquet-to-open-chain mapping hypothesis

**Status:** implemented for a pilot only; not yet validated as the physical model of the archived traces.

## Momentum-space source model

The audited FL_QOM two-step Hamiltonian is

\[
H_s(k)=
\begin{pmatrix}
M_s+J_2\cos k & J_1\sin k\\
J_1\sin k & -(M_s+J_2\cos k)
\end{pmatrix},
\qquad M_s=M-s\,\delta,
\]

with `s=0,1` for the two half-periods. The implementation uses a period
`T=2*pi`, with each half-period lasting `pi`.

## Open-chain Fourier mapping

For `L` unit cells with modes `(A_n,B_n)`, the finite single-particle matrix used
in the pilot is constructed from:

- onsite terms `+M_s` on `A_n` and `-M_s` on `B_n`;
- same-sublattice couplings `+J2/2` between neighbouring A cells;
- same-sublattice couplings `-J2/2` between neighbouring B cells;
- directed inter-sublattice couplings
  `t_AB(+1)=-i J1/2` and `t_AB(-1)=+i J1/2`, with Hermitian-conjugate B-to-A terms.

This produces the stated `sin(k)` and `cos(k)` terms in the translation-invariant
bulk convention, while open boundaries are retained in the finite chain.

## Open-system additions

The pilot adds the same dissipative and nonlinear ingredients used in the
reconstruction dynamics:

- an even cat state on the first mode and vacuum elsewhere;
- a Kerr term on the first mode;
- either low edge loss with larger bulk loss or uniform loss;
- paired disorder draws across candidate sectors;
- quantum trajectories nested within disorder realizations.

## Scientific boundary

This is a **mapping hypothesis**, not a proof that the archived N=6 static-chain
traces were generated from this driven model. It is therefore a new validation
protocol. Results from this pilot must not be merged with the archived arrays or
used to relabel them. The mapping becomes publication-grade only after:

1. the clean finite-chain Floquet spectrum and edge localization reproduce the
   selected bulk-sector labels;
2. the open-system dynamics are converged in cutoff, time step, Wigner grid and
   trajectory count;
3. matched W=1/W=0 controls are run with identical dissipation and disorder;
4. the resulting mechanism survives finite-size and loss-profile perturbations.
