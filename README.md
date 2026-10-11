# floquet-topology-synthetic-flux
We test whether a bulk invariant constrains a non-Gaussian resource under dissipation, nonlinearity, and finite size: a two-step chiral Floquet SSH chain with Kerr nonlinearity, edge-selective loss, and a non-Gaussian cat state, comparing sectors with distinct Floquet invariants under identical conditions

---

# Track A V2608 package — 2026-08-14

> **Canonical sources (single version):** `01_track_A_topo_main.tex` (main) and
> `02_track_A_topo_SI.tex` (Supplemental Material). These root-named files are the
> submission sources; the packaging files
> (cover letter, data-availability statement, suggested referees, PhySH terms,
> metadata, README_SUBMISSION, compilation commands) are kept alongside
> the canonical sources and rebuilt from them when a submission bundle is needed.

## Scientific scope

This workspace contains a bounded finite-time study of an explicitly defined
two-step chiral Floquet SSH model coupled to a finite dissipative bosonic chain. Its
primary comparison uses distinct coherent winding sectors selected from the
self-contained model audit, matched disorder realizations, and a single-site
Wigner-negativity observable. A separate micromotion-covariant pilot tests whether a
collective edge-mode population provides a more direct mechanistic observable.

The package does **not** claim:

- a topological protection theorem;
- a rigorous asymptotic decay bound;
- a universal basin or stability result;
- ultrastrong-coupling hardware validation;
- experimental feasibility of a SiN/AlN implementation;
- enhanced metrological sensitivity.

## Contents

- `01_track_A_topo_main.tex` — canonical main article;
- `02_track_A_topo_SI.tex` — canonical Supplemental Material;
- `03_track_A_topo_cover_letter.tex` — claims-aligned cover letter;
- `references.bib` — bibliography (63 entries; all cited across main+SI, verified 2026-10-11 against the rendered `.bbl` files: 63/63 cited, 0 uncited; DOI-validated, legacy metadata mismatches fixed 2026-10-02; `Plodzien2026` added 2026-10-09);
- `scripts/generate_persistence_figure.py` — deterministic reconstruction of the archived traces;
- `scripts/run_trackA_n6_independent.py` — deterministic N=6 propagation used for the independent check;
- `scripts/recompute_floquet_invariant.py` — audit of the available Floquet invariant and archived label mapping;
- `scripts/scan_floquet_phase_audit.py` — deterministic J1/J2 scan with resolution and quasienergy-gap diagnostics;
- `scripts/plot_floquet_phase_audit.py` — audit figure for the source-model scan;
- `scripts/validate_candidate_sectors.py` — resolution-ladder audit for provisional W=1/W=0 points;
- `scripts/run_matched_dynamics.py` — paired historical static-chain diagnostic with edge/uniform-loss controls;
- `scripts/run_mapped_floquet_dynamics.py` — paired dynamics for the explicit Floquet-to-open-chain mapping hypothesis;
- `scripts/audit_mapped_floquet_pilot.py` — fail-closed audit of the mapped pilot;
- `scripts/scan_open_chain_edge_modes.py` — finite-size edge-localization diagnostic;
- `scripts/validate_chiral_floquet_model.py` — independent chiral Floquet invariant and bulk-edge audit;
- `scripts/plot_chiral_floquet_audit.py` — chiral phase/edge audit figure;
- `scripts/run_chiral_floquet_dynamics.py` — dissipative dynamics for the validated chiral model;
- `scripts/run_chiral_edge_mode_pilot.py` — pre-specified edge-mode mechanism pilot;
- `scripts/analyze_edge_mode_pilot.py` — paired edge-mode/loss-current audit;
- `scripts/analyze_edge_mode_convergence.py` — Fock-cutoff convergence audit;
- `scripts/audit_chiral_floquet_pilot.py` — fail-closed dynamics-gate audit;
- `scripts/trackA_n6_2608.sbatch` — SLURM wrapper for the full N=6 run (job-name `trackA_n6_rebuild`);
- `scripts/trackA_edge_mode_pilot.sbatch` — four-condition edge-mode pilot array;
- `scripts/trackA_edge_mode_convergence.sbatch` — passive/Kerr `N_max=2,3,4` convergence array;
- `scripts/trackA_edge_mode_covariant.sbatch` — corrected micromotion-covariant pilot array;
- `scripts/run_covariant_unitary_benchmark.py` — deterministic unitary covariance and edge–bulk control;
- `scripts/trackA_covariant_unitary_benchmark.sbatch` — SLURM wrapper for the deterministic control;
- `scripts/run_passive_single_particle_benchmark.py` — deterministic loss/micromotion decomposition;
- `scripts/trackA_passive_single_particle.sbatch` — SLURM wrapper for the passive one-excitation control;
- `scripts/run_8mode_nmc500_parallel.py` — parallel N_MC=500 production runner (2/3/4 cells);
- `scripts/run_revision_scans.py` — Tier-1 revision scans (γ-ratio, multiperiod, N_max, micromotion, disorder sweep, chiral-loss);
- `scripts/analyze_revision_scans.py` — paired statistics and machine-readable summary for the revision scans;
- `scripts/run_cross_solver_benchmark.py` — `mesolve` vs `mcsolve` cross-solver validation (roadmap G5);
- `scripts/compute_bootstrap_cis.py` — bootstrap confidence intervals for production contrasts;
- `scripts/trackA_revision_4cells_nmc300.sbatch`, `scripts/trackA_revision_nmax4_3cells.sbatch`, `scripts/trackA_r1_nmax_3cells_nmc200.sbatch`, `scripts/trackA_r2_cross_solver_3cells.sbatch` — R1/R2 SLURM wrappers;
- `scripts/trackA_r4_ratio_sweep.sbatch`, `scripts/trackA_r5_chiral_loss.sbatch`, `scripts/trackA_r7_micromotion.sbatch` — R4/R5/R7 SLURM wrappers;
- `scripts/run_r3_r8_records.py` — matched-design parallel recorder for R3 (size 5) and R8 (witness rerun with terminal 3x3 blocks);
- `scripts/run_r9_longtime.py` — R9 stroboscopic long-time contrast curves (t = 0..64T);
- `scripts/predict_r3_size5.py` — freezes the out-of-sample size-5 prediction before the R3 submission (L3);
- `scripts/analyze_r3_size5.py` — scores the R3 measurement against the frozen prediction;
- `scripts/analyze_r8_witness.py` — fail-closed production-parity gate + Chabaud Fock witnesses + relative-entropy non-Gaussianity;
- `scripts/analyze_r9_longtime.py` — plateau sign persistence, stationary time and envelope-decay gap proxy;
- `scripts/verify_r3r8r9_tex.py` — fail-closed check that the R3/R8/R9 SI table numbers match their JSON sources;
- `scripts/trackA_r3_size5.sbatch`, `scripts/trackA_r8_witness.sbatch`, `scripts/trackA_r9_longtime.sbatch` — R3/R8/R9 SLURM wrappers;
- `data/trackA_persistence_summary.json` — machine-readable summary;
- `data/trackA_invariant_audit.json` — invariant audit;
- `data/trackA_heom_validation.json` — N=3 HEOM/Master Equation cross-check;
- `data/trackA_n6_independent.json` — completed deterministic N=6 finite-time check;
- `data/floquet_phase_audit.json` — completed source-model J1/J2 scan;
- `data/model_manifest.json` — explicit separation of invariant and dynamics model conventions;
- `data/candidate_sector_validation.json` — completed candidate resolution ladder;
- `data/production_run_plan.json` — paired dynamics design and staging policy;
- `data/mapped_pilot_analysis.json` — completed mapped-pilot gate audit;
- `data/open_chain_edge_audit.json` — finite-size edge-localization audit;
- `data/chiral_floquet_phase_audit.json` — independent chiral Floquet bulk/edge audit;
- `data/chiral_pilot_analysis.json` — completed chiral dynamics pilot audit;
- `data/edge_mode_unitary_benchmark.json` — machine-readable deterministic benchmark;
- `data/passive_single_particle_benchmark.json` — machine-readable passive benchmark;
- `data/edge_mode_pilot_analysis.json` and `data/edge_mode_convergence_analysis.json` — machine-readable audits;
- `data/chiral_production/chiral_floquet_{2,3,4}cells_NMC{100,300,500}.json` and `chiral_floquet_5cells_NMC300.json` — production dynamics records (2–5 cells, disorder 0–1.5, 5-cell = audit run R3);
- `data/r3_size5_prediction.json` — out-of-sample size-5 prediction, frozen before the R3 submission (L3 gate);
- `data/r3_size5_analysis.json` — R3 scoring against the frozen prediction;
- `data/r8_witness/r8_witness_{2,3,4}cells.json` and `data/r8_witness_summary.json` — audit run R8 witness rerun and certification summary;
- `data/r9_longtime/r9_longtime_{2,3}cells_NMC200.json` and `data/r9_longtime_summary.json` — audit run R9 long-time curves and summary;
- `data/chiral_production_analysis.json` and `data/chiral_production_summary.csv` — paired production analysis and summary table;
- `data/revision_scans/` and `data/revision_scans_summary.json` — R4/R5/R7 scan fragments and merged summary (γ-ratio, chiral-loss, micromotion, sweep);
- `data/cross_solver_benchmark_3cells.json` — 3-cell `mesolve`/`mcsolve` cross-validation;
- `data/convergence_checks.json` — Wigner-grid and trajectory-count sensitivity checks;
- `data/spectral_localizer_audit.json` — spectral-localizer topological certification;
- `figures/trackA_persistence_rebuilt.png` — data-driven archived-trace figure (preview only; PDF removed as unused);
- `figures/chiral_floquet_phase_audit.pdf` and `.png` — independent chiral audit figure;
- `figures/floquet_phase_audit.png` — resolution/gap audit figure (preview only; PDF removed as unused);
- `figures/trackA_n6_rebuild_summary.png` — independent N=6 check figure (preview only; PDF removed as unused);
- `figures/edge_mode_covariant_summary.pdf` and `.png` — corrected edge-mode pilot figure;
- `validation_manifest.json` — evidence and reproducibility status;
- `Track_A_data_analysis_report.md` — analysis population, provenance and limitations (canonical DAR; historical working names `DAR_TrackA.md` no longer exist);
- `PRX_QUANTUM_POPULAR_SUMMARY.txt` — nontechnical summary (PRX-series artifact;
  PRA has no popular-summary field, kept for reference);
- `CLAIM_EVIDENCE_MATRIX.md` — claim-to-evidence and claim-boundary audit;
- Historical working documents (`POSITIONING_MEMO_20261002.md`,
  `ADVERSARIAL_AUDIT_20260901.md`, `AUDIT_RESEARCH_REFINEMENTS_20261001.md`,
  and the consolidated Tracks-A–D adversarial review) were removed in the
  2026-10-09 cleanups (DAR §§26/31); their substance lives in
  `Track_A_data_analysis_report.md` §§12–13 and `CLAIM_EVIDENCE_MATRIX.md`.
  The review prompt template is
  `../PRX_PRXQuantum_adversarial_review_prompt.md`.

The revised manuscript is
now written as a bounded negative-result study rather than a computational report;
technical submission gates are tracked in `ADVERSARIAL_AUDIT_20260901.md`
(gate matrix G1–G8) and the refinement plan
`AUDIT_RESEARCH_REFINEMENTS_20261001.md`.


## Compiling the manuscript (LaTeX)

Compile **from this directory** — the source requires the bundled
`array.sty` (v2.4d) shipped here to shadow the TeX-distribution copy:
REVTeX 4.2's `tabularx` handling is incompatible with `array >= 2.5`
(``Extra \or`` errors), so a build that resolves the system `array.sty`
instead of `./array.sty` will fail with ``Newer array package detected``.
The SI preamble also aborts on that condition with an explanatory
`PackageError`, so the failure mode is explicit rather than a cryptic
cascade.

```bash
# from this directory (so ./array.sty shadows the distribution copy)
pdflatex 01_track_A_topo_main.tex
bibtex   01_track_A_topo_main
pdflatex 01_track_A_topo_main.tex
pdflatex 01_track_A_topo_main.tex

pdflatex 02_track_A_topo_SI.tex
bibtex   02_track_A_topo_SI
pdflatex 02_track_A_topo_SI.tex
pdflatex 02_track_A_topo_SI.tex
```

If you copy the sources elsewhere to build (e.g. a scratch directory or
CI), copy `array.sty` (and the `*.bib` files) alongside the `.tex` files;
compiling from a directory that contains only the `.tex`/`.bib` files is
the one configuration that reproduces the `array` incompatibility.


## Reproduction of the current figure

From this directory:

```bash
python3 -m py_compile scripts/generate_persistence_figure.py
python3 scripts/generate_persistence_figure.py
python3 scripts/generate_n6_summary_figure.py
```

The reconstruction script reads the legacy arrays in `../../code/_cache/` from the
parent `Q1_Papers` repository. Those arrays are provenance-only and are not the source
of the primary chiral-sector result. The script does not run a new QuTiP propagation
and does not fit an exponential protection law.

## Independent rerun status

The environment audit found a usable QuTiP installation at
`/home/nanaengo/miniforge3/envs/qom` (QuTiP 5.3.1, NumPy 1.26.4, SciPy 1.13.1).
The N=3 HEOM cross-check completed successfully (8/8 cases; maximum discrepancy
0.5857%), and the Floquet API audit completed. The audit also found that the two
archived nonzero labels use identical available hopping parameters and therefore
cannot currently be treated as independently verified invariant sectors. The full
N=6 deterministic rerun is submitted as SLURM job 15327; its output must be reviewed
before any manuscript claim is changed. The package is a **bounded author-review submission package** for the stated finite-time numerical result; it is not a claim of experimental or theorem-level validation.

The corrected static-chain pilot completed as SLURM job `15330`, the explicit
Floquet-mapping pilot as job `15334`, and the independent chiral-model pilot as job
`15336`. The chiral bulk/edge audit passed, but the dissipative Wigner-persistence
gate did not yet pass at 6 modes or under uniform loss. The original large production
run was therefore blocked rather than scaled blindly.

## Chiral production checkpoint — job 15338

A bounded confirmatory production was subsequently executed with the independently
audited chiral Floquet SSH model (SLURM array job `15338`). The primary production inventory is now updated to the complete JSON outputs used by the DAR:

- 4 modes (2 cells): `N_MC=500`, 8,000 records, 4,000 paired units;
- 6 modes (3 cells): `N_MC=500`, 8,000 records, 4,000 paired units;
- 8 modes (4 cells): `N_MC=300`, 4,800 records, 2,400 paired units; the earlier `N_MC=100` file is retained as a nested historical comparison.

All outputs passed structural auditing (sector pairing, finite values, declared
seeds and record counts). The combined analysis is in
`data/chiral_production_analysis.json`, with the tabular summary in
`data/chiral_production_summary.csv` and the figure in
`figures/chiral_production_contrast.{png,pdf}`.

The edge-loss contrast was positive for 2 cells (approximately 0.013–0.017 across
the tested disorder levels), but negative for 3 cells (approximately −0.0020 to
−0.0023) and also negative for 4 cells (approximately −0.0018 to
−0.0026; the current `N_MC=300` superset gives the primary interval). Uniform-loss contrasts were near zero
or slightly negative at all sizes. Thus, the production result does **not** support
a size-robust claim that the certified Floquet sector protects local Wigner
negativity. It is retained as a bounded negative/size-sensitivity result, not as a
protection theorem, metrological claim, or experimental prediction.

The output is descriptive finite-time evidence only. The original labels, the
assumed reduced-model parameters, the local Fock cutoff, and the finite disorder
ensemble remain part of the declared scope. No manuscript claim is upgraded beyond
this boundary.

## Next implementation: edge-mode mechanism pilot

The production result showed that local site-0 Wigner negativity is not size-robust.
The next authorized computation is therefore the pre-specified edge-mode pilot, not
another large Monte Carlo run. It uses the same chiral Floquet backbone, the
collective operator `b_edge`, paired loss profiles, common edge-site and
mode-resolved cat preparations, and passive/Kerr controls.

The pilot was executed as SLURM job `15343` and the cutoff panel as job `15344`.
Both arrays completed with empty error logs and structurally complete outputs.The four pilot files contain 160 records each; the six convergence files contain six
records each. The paired audits are `data/edge_mode_pilot_analysis.json` and `data/edge_mode_convergence_analysis.json`.

The structural gate passed, but the passive-reference residuals were 0.03883,
0.02821 and 0.03557 for `N_max=2,3,4`, respectively. They were not monotone and no
declared numerical-error budget was met. The passive-bound, cutoff-convergence and
topology-by-loss mechanism gates therefore remain **NOT PASSED**. The collective-mode
Wigner transform is still not validated; local site negativity is not its substitute.
All pilot parameters remain normalized, uncalibrated and unoptimized. No larger run
or manuscript protection claim is authorized from this result.

## Corrected micromotion-covariant pilot — job 15353

The fixed-frame ambiguity was corrected before another production: the standard
number operator was implemented and the edge vector was transported through the
coherent Floquet micromotion. A local smoke test passed before submission. Array job
`15353` then completed all four conditions with empty stderr logs and 160 records per
condition. The audit is `data/edge_mode_covariant_analysis.json`.

At zero disorder and edge-low loss, the common-site passive contrast was only
0.000388 (95% interval [0.000358, 0.000418]); the topological edge-mode preparation
showed a small negative contrast of −1.59×10⁻⁵ (95% interval [−1.82×10⁻⁵,
−1.36×10⁻⁵]). The corrected passive residuals reached 0.03956 for the common-site
preparation and 0.38257 for the mode-resolved preparation. The passive bound and
mechanism gates therefore remain **NOT PASSED**. This correction is a diagnostic
result, not evidence of protection. The subsequent deterministic unitary
micromotion/edge–bulk benchmark passed for 4, 6, 8 and 12 modes, confirming the
observable/frame correspondence. A one-excitation passive benchmark is exact under
uniform loss but shows a small structured-loss residual and complementary-space
population. It explains the open-system reference failure at model level, but does
not validate protection; additional Monte Carlo scaling remains unjustified.

## Reviewer-response updates (2026-08-14)

- **Eight-mode result reframed**: the zero-disorder contrast at eight modes is
  −0.00141 with a bootstrap 95% interval [−0.00392, +0.00097] (n=100) that crosses
  zero; abstract, Results, and cover letter now describe attenuation toward zero
  rather than a negative value. Bootstrap percentile intervals (10⁴ resamples over
  disorder realizations) replace the normal-approximation display; contrasts are
  descriptive and not p-value tested.
- **Introduction**: explicit novelty statement and operational definition of
  protection (three prespecified criteria, including a 10⁻³ passive-reference
  budget, which both preparations violate by ~40× and ~380×).
- **SI Sec. S9 (new)**: convergence diagnostics — local-contrast cutoff dependence
  (+0.0154/+0.0109/+0.0249 at N_max = 2/3/4, non-monotone, sign robust at two
  cells), Wigner-grid integral within ±3% of the finest grid, trajectory-count
  sensitivity ~4% on a fixed realization. Provenance:
  `data/convergence_checks.json` (31 s, 240 records),
  `scripts/run_convergence_checks.py`.
- **siunitx/cleveref hygiene** applied to main, SI, and cover letter.

## V2608 rename and roadmap implementation (2026-08-14, second pass)

- **Rename**: package renamed from `track_A_topo_rebuilt_20260814` to
  `track_A_topo_2608`; canonical sources were named `01_track_A_topo_main_2608.tex` and
  `02_track_A_topo_SI_2608.tex` at that point (the package was later re-based
  onto the root names `01_track_A_topo_main.tex` / `02_track_A_topo_SI.tex`,
  which are the current canonical sources); historical audit files consolidated into
  `../PRX_PRXQuantum_adversarial_review_tracks_A_B_D_20260822.md`; all provenance
  paths in data JSONs and sbatch scripts updated; obsolete artifacts
  removed (stray SLURM err/out, `__pycache__`, legacy figures
  `floquet_phase_audit`, `trackA_n6_rebuild_summary`, `trackA_persistence_rebuilt`
  — all regenerable from retained data + scripts).
- **Figure 1 model schematic (roadmap G11)**: new `figures/model_schematic.pdf/png`
  (chain + drive protocol + loss channels + even-cat observable), generated by
  `scripts/plot_model_schematic.py`, referenced as `\cref{fig:schematic}` in the
  main text. Main is 15 pages (5 figures).
- **Cross-solver benchmark (roadmap G5, item 6)**: `scripts/run_cross_solver_benchmark.py`
  integrates the exact Lindblad master equation (QuTiP `mesolve`) alongside the
  Monte Carlo solver (`mcsolve`) on the two-cell subset — 8 conditions (2 sectors ×
  2 loss profiles × disorder 0.0/0.5), ntraj=200. Mean Uhlmann fidelity
  0.99957, max trace distance 1.93e-2, max occupation RMSE 6.7e-3. Results in
  `data/cross_solver_benchmark.json`; new SI S9 subsection. SM is 19 pages
  (§S1–S11, Figs. S1–S4).
- **Collective-mode Wigner transform (roadmap G5 extension)**: validated bosonic
  mode rotation (`scripts/run_collective_wigner.py`) closes the collective-mode
  gate with a measured null result — collective edge-mode Wigner negativity is
  not resolved above the tested numerical floor in all 80 matched records (edge-mode Wigner positive
  semidefinite), while the local site-0 negativity is positive under
  edge-selective loss (R14; SI §S9, Fig. S2).
- **Run manifest (roadmap G1)**: `data/run_manifest.json` documents seeds,
  commands, jobs, and claim boundary for all 14 runs (incl. R13
  spectral-localizer certification and R14 collective-mode Wigner).
- **Environment lock (roadmap G12)**: `data/qom_environment_lock.txt` (conda
  export, 212 packages) and `data/qom_pip_lock.txt` (pip freeze, 59 packages).

## R3/R8/R9 integration, references, and positioning pass (2026-10-02)

- **R3/R8/R9 integrated** into main + SI (size-5 out-of-sample table, grid-free
  witness table, 64-period table) with `scripts/verify_r3r8r9_tex.py` as a
  fail-closed table↔JSON gate (currently FAILS: 0).
- **Figure deduplication**: the two SI copies of main Figs. 4 and 6 were removed;
  the SI now carries exactly two unique figures (S1 protection bound,
  S2 spectral-localizer profiles) and refers to the main-text figures in prose.
- **References expanded 28 → 41**, all cited across main+SI, metadata verified
  against CrossRef/OpenAlex; two metadata bugs fixed (`Diehl2011` DOI
  nphys2081 → nphys2106; `Bardyn2013` title restored to "Topology by
  Dissipation"); three uncited dead entries removed; `validate_citations.py`
  reports 0 errors / 0 duplicates.
- **Positioning moves M1–M5 implemented** per `POSITIONING_MEMO_20261002.md`:
  reconciliation paragraph on prior non-Gaussian/topological protection claims,
  invariant-null clause in the abstract, paired non-Gaussianity contrast
  (+0.1316±0.0025) surfaced in Results, R5 strong-symmetry numbers in the
  Discussion, and the cover letter rewritten to the composite identity
  ("real, certified, and predictable, but not protected"; stale 1.8× power
  figure corrected to 1.35×).
- **Package cleanup**: stale draft `04_popular_summary.txt` (unreferenced,
  "two to eight cells") removed in favour of the canonical
  `PRX_QUANTUM_POPULAR_SUMMARY.txt` (enriched with the out-of-sample prediction
  and grid-free certification; 173/250 words); `README_SUBMISSION.md` figure
  inventory corrected to Figs. 1–7 + SI S1–S2.
- **Gates**: main 21 pp / SI 37 pp / cover 2 pp; 0 LaTeX errors, 0 undefined
  refs, 0 `??`, overfull unchanged (2 pre-existing), forbidden-phrase sweep 0
  hits across main+SI+cover.
- **M6 figure economy implemented**: `fig:phase` demoted from the main text to
  SI Fig. S1 (main now 6 figures, SI 3; main §II points to `SI-fig:phase`),
  with `README_SUBMISSION.md` figure inventory updated accordingly.

### 2026-10-02 (desk-reject audit + SI overhaul)

- **Skills applied from `/home/nanaengo/skills`**: `peer-review`,
  `venue-templates`, `scientific-writing`, `citation-management`. PRX Quantum
  author instructions fetched 2026-10-02 (journals.aps.org/prxquantum/authors):
  research articles have no length limit; all references cited in the
  Supplemental Material must be listed in the main reference list; reference
  titles encouraged; popular summary ≈150 words.
- **Desk-reject audit fixes (main)**: APS SM-reference violation repaired
  (`Johansson2013` was cited only in the SI — the solver sentence in §III now
  carries it); abstract compressed 295 → 271 raw words while keeping every
  verified number; author-contributions placeholder replaced with the
  per-author statement supplied by the team (JD/PD simulations and data
  analysis; AHAA/AE/NM theory and feedback; SGNE conception, supervision, and
  direction; PD/SGNE manuscript writing; all authors approved).
- **SI data-consistency repair**: stale pre-upgrade counts in §S4/§S5,
  `tab:paired-comparison`, and `tab:parameters` (`N_MC=300` / 2400 units /
  range −0.0026…−0.0018 for the eight-mode production) corrected to the
  production truth `N_MC=500` / 4000 units / −0.0023…−0.0020, verified against
  `data/chiral_production_summary.csv` and `data/power_analysis.json` (the
  2026-09-01 audit's 300→500 upgrade is now reflected everywhere).
- **"fifth size" → "five-cell size"** sweep across main abstract/§III/
  Discussion/Conclusion, SI §S9 title+body, `tab:scope-limitations`, and the
  popular summary — removes the off-by-one ordinal (four sizes are listed).
- **SI references 10 → 15 keys**: added `Dalibard1992` (MCWF method),
  `Efron1979` (bootstrap), `Wilson1927` (Wilson intervals) with Crossref-
  verified DOIs; the raw `arXiv:2410.24176` placeholder in the localizer bound
  replaced by `\cite{WongCerjanCole2024Quantitative}`; `Xuereb2014` added to
  the hardware mapping; bib 41 → 44 entries, all 15 SI keys also cited in the
  main text (APS SM rule), `validate_citations.py`: 0 errors / 0 duplicates
  (1 recommended-field warning: Efron page range absent from both Crossref and
  OpenAlex, left unfilled by no-fabrication policy).
- **SI prose**: scope section expanded into a navigation paragraph (full size
  set + extended analyses), hardware section's sentence fragment fixed.
- **Measurement-chain paragraph (SI §S12, answering the experimental-realism
  audit)**: site Wigner reconstruction via parity-based direct measurement
  (`Lutterbach1997`), displacement-and-readout primitive demonstrated with
  111-photon cat states in circuit QED (`Vlastakis2013`), reservoir-stabilized
  manifolds for the engineered bulk loss (`Touzard2018`), and homodyne
  quadrature readout on the optomechanical side (`Aspelmeyer2014`); main-text
  Discussion now cites all four (APS SM rule). Bib 44 → 48 entries,
  `validate_citations.py`: 0 errors / 0 duplicates.
- **HPC convergence runs (SLURM, partition `production`)** — **complete**:
  job **15914** `scripts/trackA_nmax3_production.sbatch` — cells {2,3,4} at
  `N_max=3`, `N_MC=500`, production seeds (matched design; output
  `data/chiral_production_nmax3/`) — finished the 2-/3-cell files 2026-10-02;
  the 4-cell phase was infeasible in full space (dim 65536, ~500 h needed) and
  was superseded 2026-10-03 by job **15916** after the runner gained an exact
  total-quanta Q<=2 subspace reduction (dim 45; invariance + record-level
  validation in `scripts/validate_subspace_reduction.py`), which completed all
  8000 records in 49 s. Job **15915** `scripts/trackA_nmc2000_8mode.sbatch` —
  8 modes at `N_MC=2000`
  (4× the production sample, nested samples 0–499; output
  `data/chiral_production/chiral_floquet_4cells_NMC2000.json`).
  `run_8mode_nmc500_parallel.py` has row-level checkpointing
  (`.partial.jsonl` resume, canonical record order preserved — physics and
  seeds untouched) and skips finalized outputs.
  **Analysis done 2026-10-03**: `scripts/compare_nmax_truncation.py` →
  `data/nmax2_nmax3_paired_comparison.json` (paired N_max=2 vs 3 shift:
  excludes zero at 2 cells, −0.00540 [−0.00615, −0.00467]; consistent with
  zero at 3/4 cells) and `scripts/compute_bootstrap_cis.py --prod-dir
  data/chiral_production_nmax3 --out data/bootstrap_contrast_cis_nmax3.json`
  (three-size N_max=3 CIs, sign pattern confirmed at N_MC=500); DAR §11
  updated first, then SI `sec:nmaxpanel`/`tab:nmaxpanel` and main-text
  limitations narrowed accordingly.
- **Gates**: all three PDFs rc=0; main 21 pp / SI 40 pp / cover 2 pp; 0 LaTeX
  errors, 0 undefined references, 0 `??`; overfull unchanged (2 pre-existing);
  `verify_r3r8r9_tex.py` FAILS: 0; forbidden-phrase sweep 0 hits.

## Journal retarget — Physical Review A (2026-10-03)

- **Decision:** Track A target switched from PRX Quantum to **Physical Review A**
  (author decision: Q1 in Atomic/Molecular/Optics, APC-free via the APS
  subscription route; PRX Quantum is gold OA, 2026 APC $3590).
- Sources switched to the REVTeX **`pra`** class option (main + SI); gates
  re-run after retarget: main 21 pp / SI 40 pp / cover 1 p; 0 errors,
  0 undefined refs, 0 `??`; overfull baseline unchanged (single page-1
  REVTeX vbox); `verify_r3r8r9_tex.py` FAILS: 0.
- Retargeted documents: cover letter (addressee + article-type wording),
  `metadata.md`, `README_SUBMISSION.md` (journal/article-type checkboxes,
  popular summary marked not-required for PRA, hybrid APC note),
  `05`/`06`/`07` headers, `compilation_commands.txt`
  (`prx` → `pra` class option), `validation_manifest.json` (checklist pointer
  fixed — the previous target file did not exist; `target_journal` recorded).
- Historical audit/positioning documents (`ADVERSARIAL_AUDIT_*`,
  `POSITIONING_MEMO_*`, `AUDIT_*`, literature audits) keep their original
  PRX Quantum framing as dated records.

## Positioning re-interrogation and mitigations (2026-10-03, second pass)

- **Trigger:** re-interrogate positioning after a deep web search with the
  appropriate skill; mitigate meticulously (report-before-tex).
- **Skill/fallback provenance:** `research-lookup` skill loaded; its parallel
  backends are unavailable (`parallel-cli` absent, `PARALLEL_API_KEY` absent,
  OpenRouter Perplexity fallback HTTP 402 — evidence
  `data/positioning_reinterrogation_20261003/skill_fallback.log`,
  context in same directory). Documented fallback: **10 built-in web-search
  passes** (2026-10-03), all hits metadata-checked before use.
- **Findings:** no scoop and no contradiction — no prior/concurrent study
  reports local Wigner-negativity / non-Gaussian-resource tracking (or its
  failure) under a bulk Floquet invariant in an open chain, and no
  topology→non-Gaussian-resource inference-limit paper exists. PRA fit
  re-confirmed (PRA publishes WN witnessing, non-Gaussian resource no-gos,
  driven-dissipative bosonic chains). Full pass table and mitigation register:
  `DAR §12`.
- **Mitigations executed (exactly two `.tex` edits + one bib entry):**
  M8 — priority hedge "To our knowledge," before the intro's absence claim;
  M9 — cite `Raza2026` (Raza et al., *Phys. Rev. A* **114**, 013708 (2026),
  CrossRef-verified) in the Discussion parallels paragraph: Kerr+loss
  Wigner-negativity generation with no topological structure, which
  *decouples* the resource from topology. Five further citation additions
  were considered and skipped with rationale (DAR §12.4, M10–M14).
- **Gates (2026-10-03 positioning pass):** SI→main rebuilt; 0 `!` errors,
  0 undefined references (SI log's `OMS/cmtt/m/n undefined` is a pre-existing
  font-shape warning, not a reference warning), 0 `??` in all three PDFs;
  overfull baseline unchanged (1, page-1 REVTeX vbox); `verify_r3r8r9_tex.py`
  FAILS: 0; forbidden-phrase sweep 0; AI-tell word sweep 0;
  `validate_citations.py` 0 errors / 0 duplicates (only the pre-existing
  Efron1979 `pages` warning and cross-file unused warnings);
  **main 22 pp** (was 21 — the added reference pushed the bibliography onto
  one extra page; PRA has no page cap) / SI 40 pp / cover 1 p;
  M8/M9 confirmed rendered via `pdftotext -raw`.
- Checksums regenerated and verified (**478** entries: the previous 477 plus
  the new `data/positioning_reinterrogation_20261003/context.json`; the skill
  failure log is excluded by the standing `*.log` exclusion).

## Three-readings pass (2026-10-03: authors / reviewers / editor)

- **Trigger:** user instruction for a careful triple reading of the
  submission. Skill: `peer-review` (local-only working review of our own
  manuscript; no external review service). Full findings and decisions:
  `DAR` §13 (written before any `.tex` edit,
  per report-before-tex).
- **Author reading — 10 findings (A1–A10), all fixed:**
  A1 — main claimed χ/ω_d≈0.08 while the production JSON gives
  `kerr=0.5`, `period=2π` ⇒ χ/ω_d=0.5 (the 0.08 was 0.5/2π);
  A2 — SI hardware section's stated rule (Γ→Γ·ω_d) contradicted its own
  arithmetic (×ω_d/2π) and its "rates in units of 1/T" line; the mapping
  was recomputed consistently ×ω_d on both platforms (verified: this is
  the only convention that reproduces the simulated per-period decay and
  the section's own Q_edge=10⁴); ratios χ/J_B, γ_bulk/χ, 5×10⁻⁴ unchanged,
  no claim scope change (section is explicitly "not a calibration");
  A3/A4 — two stale hard-coded main-figure numbers in the SI (Fig. 6/4)
  converted to `\cref{MAIN-fig:covariant}` / `\cref{MAIN-fig:collective-wigner}`
  (actual numbers: Fig. 5 / Fig. 3);
  A5 — `tab:boundary` "contrast vanishes" contradicted the text's own
  resolved small negative (−0.00081±0.00005) → "collapses to a small
  negative value with no sign change across sizes";
  A6 — `fig:production` caption (b) scoped to the body's exact claim
  (zero for d≤0.5, small offsets at stronger disorder);
  A7 — "above any plausible numerical tolerance" → the prespecified
  10⁻³ residual threshold with the 40–380× factors (main L224 and L260);
  A8 — intro now names all sizes (four to six, eight, and ten modes);
  A9 — "changes sign at three and four cells" → "negative at three,
  four, and five cells"; A10 — remaining Unicode dashes → TeX `---`/`--`.
- **Reviewer reading (RV1 action):** abstract now states the per-stratum
  counts "(500/500/2000 per disorder level and loss profile)" behind the
  "at least 4000 matched disorder units" claim. All other referee checks
  (cutoff caveat, power, extreme loss ratio, operational protection
  definition, statistical language) were already carried — see §13.2.
- **Editor reading:** venue fit, abstract/title/figure structure,
  declarations, cover letter all pass; the two stale SI figure numbers
  (ED5) were the only defects. See §13.3.
- **Gates (2026-10-03 three-readings pass):** SI→main rebuilt; 0 `!`
  errors, 0 undefined references, 0 `??` in both PDFs; overfull baseline
  unchanged (1 unique, the page-1 REVTeX vbox — the log counts it twice,
  once per pass); `verify_r3r8r9_tex.py` FAILS: 0; forbidden-phrase sweep
  0; AI-tell sweep 0; `validate_citations.py` main 0 errors / 0
  duplicates (SI: 1 error, the pre-existing `@ifpackagelater` preamble
  false-positive verified present at git HEAD; Efron1979 `pages` warning
  unchanged); **main 22 pp / SI 40 pp / cover 1 p** (unchanged);
  main text 7122 words; every edited claim re-verified in the rendered
  PDFs via `pdftotext`.
- Checksums regenerated and verified (**478** entries, 0 failures;
  the `logs_build_*.log` files are excluded by the standing `*.log`
  exclusion).
- **Provenance follow-up (DAR §13.5):** the two numbers previously flagged
  unverified were traced to source — main L89 gaps (2.1991 / 0.3142) =
  `minimum_zero_gap` / `minimum_pi_gap` in
  `data/chiral_floquet_phase_audit.json`; main L224 covariant triplets =
  mean [95% interval] in `data/edge_mode_covariant_analysis.json`
  (common-site `files/0/summaries/0`, mode-resolved `files/2/summaries/0`).
  Editorial fix A11: both triplets now use the established
  `mean [low, high]` bracket convention with an explicit
  "(mean and 95% interval)" label. Rebuilt: 0 `!` / 0 undefined / 0 `??`,
  overfull baseline 1, `verify_r3r8r9_tex.py` FAILS: 0, main 22 pp.
  Checksums regenerated after this entry: **478 entries, 0 failures**
  (run `sha256sum -c checksums.sha256` after any further edit).

## Five-cell ntraj200 expansion (2026-10-10/11)

- **Run**: `data/solver_contrast_check_5cells_n200_nmc100.json` (r3 relaunch,
  n=100 matched draws, ntraj=200, edge d=0; log
  `logs/contrast_5cells_n200_nmc100_r3.log`, 200/200) — completed
  2026-10-10: contrast $-0.000288\pm0.000151$ (SE), 95% CI
  $[-0.000584,+0.000008]$, crossing zero (z=-1.91): the pre-committed
  honest unresolved null (DAR §34).
- **Ladder refit** (`scripts/analyze_solver_contrast_ladder.py`, dedupe
  keeps max-n): five-cell model B $C=-0.000393$, model A
  $C=-0.000183\pm0.000068$; SI `tab:trajladder` + prose, main Results
  range, and matrix C25/C29 refreshed accordingly. Five cells remain
  solver-unresolved, resting on validated extrapolation.
- **Gates**: main 26 pp / SI 49 pp / cover 1 p (rebuild 0 errors,
  0 undefined); `verify_r3r8r9_tex.py` FAILS: 0; `checksums.sha256`
  541/541 OK (new JSON appended).
