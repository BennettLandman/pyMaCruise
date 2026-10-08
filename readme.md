# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-08

The restartable synthetic base workflow completes all seven stages on the brain-resolution grid, including two naturally converged TOADS passes, CRUISE core, MaSeg refinement, and scientific surfaces/thickness. The manifest passes all stage checks; all 31 stage NIfTIs are finite RAS. Same-settings resume skips all seven stages. Peak sampled RSS was 8.433 GB for the run and 0.124 GB for resume. A current-run selected Java TGDM comparison passes all three agreed surface gates. The production surface extractor now matches JIST's connectivity-18 case templates; on the saved base level sets all three meshes have exact ordered triangle correspondence and a maximum paired-vertex displacement of 7.2×10⁻⁶ mm. See the [base workflow report](docs/reports/base-workflow-seven-stage-brain-resolution-2026-10-08.md), [current-run TGDM comparison](docs/reports/base-workflow-current-run-java-tgdm-2026-10-08.md), and [JIST template comparison](docs/reports/jist-connectivity18-case-template-port-2026-10-08.md). The naturally converged 64×80×64 synthetic plus workflow completes all nine stages, with LesionTOADS converging at iteration 6; its peak sampled RSS was 0.824 GB, and resume peaked at 0.236 GB. All three saved plus level sets also match ordered JIST faces exactly with maximum paired-vertex displacement below 0.000001 mm. These are synthetic integration and selected-source results, not whole-workflow MATLAB/JIST equivalence or real-MRI validation.

An isolated installed wheel also completed all seven base stages from outside the source checkout, loaded bundled atlas assets, and passed all stage checks. Independent inspection verified 46 output hashes, 31 finite RAS NIfTIs, and four meshes; resume skipped every stage. TOADS was intentionally capped at two iterations for this packaging check, so convergence evidence remains the separate source-checkout run. See the [installed-wheel report](docs/reports/installed-wheel-seven-stage-base-workflow-capped-2026-10-08.md).

The isolated wheel also completed the nine-stage synthetic plus workflow and naturally converged at LesionTOADS iteration 6; an independent audit verified 86 manifest-referenced file hashes, 62 finite RAS NIfTIs, three meshes, and a successful resume. Both wheel runs are synthetic packaging evidence, not whole-workflow Java/MATLAB equivalence or real-MRI validation. See the [plus installed-wheel report](docs/reports/installed-wheel-plus-workflow-2026-10-08.md).

P18ck is reopened: its earlier pyramid comparison exercised a helper transcription, not the actual VandyCRUISE pyramid method. The new direct source harness matches final transform and lambda exactly; per-level energy/stop-state equivalence remains under investigation. See the [actual-source report](docs/reports/toads-registration-pyramid-actual-source-2026-10-08.md).

The full suite on tested tree `1bb0ab8` passed **1,500 tests, 0 failures, 0 skips, and 11 expected RAS-orientation warnings** in **1,815.06s**, peaking at **13.010 GB RSS** under a 16 GB watchdog; see the [suite report](docs/reports/full-pytest-2026-10-08-jist-connectivity18.md). Remaining work includes source-composite comparisons and refinements, plus cross-grid surface-position differences, paper-based code annotations, and real-MRI validation. The available human T1 still needs its matching user-provided brain mask.
## Scope and target

- **Base workflow:** T1-weighted MRI, MaSeg labels, and a user-provided brain mask on the T1 grid.
- **Plus workflow:** T1, MaSeg, matching brain mask, and FLAIR, following the selected `run_maCRUISEplus_v1` profile in [docs/plus-baseline.md](docs/plus-baseline.md).
- **Reference source:** MATLAB MaCRUISE V3 and the supplied VandyCRUISE algorithms. Numbered historical runners are bug-fix versions; versions before V3 contain known bugs.
- **Goal:** run the selected workflows without MATLAB or a MIPAV/JIST JVM runtime. ANTsPy/ANTs is the selected external backend for ANTs operations.
- **Excluded:** visualization porting and FreeSurfer requirements. A post-port discussion will select modern Python visualization tools.

## Implementation principles

1. Use existing Python libraries for general image I/O, array arithmetic, registration, filtering, morphology, and mesh operations when their behavior fits the source; port source-specific behavior with direct evidence where it does not.
2. Preserve useful MaCRUISE outputs, labels, data types, transforms, and source parameters. Treat image calculation as a generic geometry-aware N-dimensional array operator.
3. Save image inputs, fixtures, intermediates, and outputs as compressed NIfTI (`.nii.gz`). Report spatial outputs in RAS; warn and record reorientation of non-RAS inputs.
4. Build and validate on deterministic synthetic data before using the supplied human scan. Keep synthetic fixtures and generators in `tests/data/synthetic/`.
5. Report source-equivalence evidence, Python integration regressions, and real-data validation as distinct evidence categories. Do not infer whole-workflow equivalence from component matches.

The runtime stack is NiBabel, NumPy, and SciPy; pytest, build, and scikit-image are development dependencies, and `antspyx` is optional for ANTs-backed operations. The license/attribution inventory and follow-up TTO work are in [license_info.md](license_info.md). All atlas images needed at runtime are bundled in this repository as regular Git files.

## Inputs and outputs

The command-line input, assumption, and output contract is documented in [docs/cli-help.md](docs/cli-help.md). Coordinate and NIfTI rules are in [docs/image-coordinate-conventions.md](docs/image-coordinate-conventions.md). Numerical acceptance thresholds are in [docs/validation-tolerances.md](docs/validation-tolerances.md).
