# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-04

The restartable brain-resolution synthetic base workflow has now run through all six stages, and a same-parameter resume verified and skipped every stage. Outputs are finite RAS NIfTI on the native grid; peak sampled RSS was 8.308 GB under a 16 GB watchdog. The first TOADS pass hit its 50-iteration cap without reaching the energy threshold, while the second converged at iteration 20. The manifest therefore records `status=complete`, `converged=false`, and `source_equivalence_claim=false`. This is full-route synthetic integration evidence, not whole-workflow source equivalence or clinical validation. Details and reproduction commands are in the [brain-resolution base-run report](docs/reports/base-workflow-brain-resolution-converged-2026-10-04.md).

The selected Java multigrid GVF comparison passes at full crop size (maximum component error `6.77e-6` against a `5e-5` gate). The latest full-suite result remains **1,141 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 801.80s** on `ea31121`; no pytest suite was run for the new runtime integration. Remaining base-workflow evidence includes resolving the first-pass TOADS convergence/source mismatch and stage-wise Java/MATLAB output comparisons. The plus workflow, real-MRI validation, and published-paper documentation/annotations remain open. See [plan.md](plan.md) for the detailed steps and current blockers.
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
