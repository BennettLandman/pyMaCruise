# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-08

The restartable synthetic base workflow completes all seven stages on the brain-resolution grid, including two naturally converged TOADS passes, CRUISE core, MaSeg refinement, and scientific surfaces/thickness. The manifest passes all stage checks; all 31 stage NIfTIs are finite RAS. Same-settings resume skips all seven stages. Peak sampled RSS was 8.433 GB for the run and 0.124 GB for resume. A current-run selected Java TGDM comparison passes all three agreed surface gates; direct JIST surface extraction comparison passes sampled position limits but confirms the Lewiner connectivity deviation. See the [base workflow report](docs/reports/base-workflow-seven-stage-brain-resolution-2026-10-08.md), [current-run TGDM comparison](docs/reports/base-workflow-current-run-java-tgdm-2026-10-08.md), and [current-run JIST surface comparison](docs/reports/base-workflow-current-run-jist-surfaces-2026-10-08.md). The naturally converged 64×80×64 synthetic plus workflow also completes all nine stages on the current code fingerprint, with LesionTOADS converging at iteration 6; its peak sampled RSS was 0.824 GB, and resume peaked at 0.236 GB. The current plus level sets have matching JIST/Python component counts and Euler characteristics with sub-0.0001 mm maximum vertex-to-triangle distances, though face connectivity differs. See the [plus workflow report](docs/reports/plus-workflow-current-code-64-grid-2026-10-08.md) and [current-run surface comparison](docs/reports/plus-workflow-current-run-jist-surfaces-2026-10-08.md). These are synthetic integration and selected-source results, not whole-workflow MATLAB/JIST equivalence or real-MRI validation.

The full suite on code-equivalent commit `5a412b7` passed **1,243 tests, 0 failures, 0 skips, and 11 expected RAS-orientation warnings** in **1,743.86s**, peaking at **14.961 GB RSS**; see the [suite report](docs/reports/full-pytest-2026-10-08-base-surface-thickness.md). Remaining work includes source-composite comparisons and refinements, plus cross-grid surface-position differences, paper-based code annotations, and real-MRI validation. The available human T1 still needs its matching user-provided brain mask.
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

The core stack is NiBabel, NumPy, SciPy, scikit-image, and pytest; `antspyx` is optional and required for ANTs-backed operations. The license/attribution inventory and follow-up TTO work are in [license_info.md](license_info.md). All atlas images needed at runtime are bundled in this repository as regular Git files.

## Inputs and outputs

The command-line input, assumption, and output contract is documented in [docs/cli-help.md](docs/cli-help.md). Coordinate and NIfTI rules are in [docs/image-coordinate-conventions.md](docs/image-coordinate-conventions.md). Numerical acceptance thresholds are in [docs/validation-tolerances.md](docs/validation-tolerances.md).
