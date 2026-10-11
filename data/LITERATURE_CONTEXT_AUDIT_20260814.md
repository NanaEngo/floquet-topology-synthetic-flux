# Literature-context audit — Track A

**Date:** 2026-08-14  
**Purpose:** document the literature checks used to refine the Introduction and Discussion. These references provide context and methodological boundaries; none is evidence for the numerical results in this package.

## Verified references used in the manuscript

| Key | Verified bibliographic record | Role in the manuscript | Evidence source |
|---|---|---|---|
| `Asboth2014` | J. K. Asbóth, B. Tarasinski, and P. Delplace, *Chiral symmetry and bulk-boundary correspondence in periodically driven one-dimensional systems*, **Phys. Rev. B 90**, 125143 (2014), DOI `10.1103/PhysRevB.90.125143`. | Symmetric-time-frame winding construction and coherent chiral Floquet sectors. | APS DOI page inspected during the literature search. |
| `Mori2023` | T. Mori, *Floquet States in Open Quantum Systems*, **Annual Review of Condensed Matter Physics 14**, 35–56 (2023), DOI `10.1146/annurev-conmatphys-040721-015537`. | Periodically driven open-system master equations and the distinction between coherent Floquet theory and dissipative dynamics. | Annual Reviews record identified by search; arXiv record `2203.16358` confirms title, author, journal DOI and scope. |
| `Xu2022micromotion` | P. Xu, W. Zheng, and H. Zhai, *Topological Micromotion of Floquet Quantum Systems*, **Phys. Rev. B 105**, 045139 (2022), DOI `10.1103/PhysRevB.105.045139`. | Motivation for treating micromotion as part of the observable construction. | arXiv record `2106.14628` confirms journal reference and DOI. |
| `Paszko2024` | D. Paszko, D. C. Rose, M. H. Szymanska, and A. Pal, *Edge Modes and Symmetry-Protected Topological States in Open Quantum Systems*, **PRX Quantum 5**, 030304 (2024), DOI `10.1103/PRXQuantum.5.030304`. | Boundary condition that open-system edge protection depends on the symmetry and structure of the dissipators. | arXiv record `2310.09406` confirms authors, title, journal reference, DOI and abstract scope. |
| `Chabaud2021witnessingwigner` | U. Chabaud, P.-E. Emeriau, and F. Grosshans, *Witnessing Wigner Negativity*, **Quantum 5**, 471 (2021), DOI `10.22331/q-2021-06-08-471`. | Distinguishes a finite numerical Wigner integral from operational witnesses/certification. | Quantum article page inspected; page supplies BibTeX, DOI and abstract. |
| `Zaw2024` | L. H. Zaw, *Certifiable Lower Bounds of Wigner Negativity Volume and Non-Gaussian Entanglement with Conditional Displacement Gates*, **Phys. Rev. Lett. 133**, 050201 (2024), DOI `10.1103/PhysRevLett.133.050201`. | Context for certifiable lower bounds, without claiming that the present calculation implements that protocol. | APS search record identified by exact title and DOI. |
| `Sieberer2025` | L. M. Sieberer, M. Buchhold, J. Marino, and S. Diehl, *Universality in Driven Open Quantum Matter*, **Rev. Mod. Phys. 97**, 025004 (2025), DOI `10.1103/RevModPhys.97.025004`. | Recent context for driven-dissipative matter and why open-system observables need not follow the coherent backbone alone. | APS search record and arXiv record `2312.03073` identified by exact title, authors, journal and DOI. |

## Claim boundary

The manuscript uses these sources only for methodological context. It does **not** claim that:

- the present finite chain realizes any of the regimes established in the cited open-system studies;
- the present local Wigner integral is an experimentally certified witness;
- the cited literature validates the present edge-mode passive reference;
- the cited literature upgrades the present result to topological protection.

The numerical claims remain bound to the local production files, deterministic controls,
and claim–evidence matrix in the package root.
