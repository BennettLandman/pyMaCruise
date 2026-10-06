# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-06

The restartable brain-resolution synthetic base workflow completes all six stages, converges both TOADS passes (iterations 39 and 25), writes 30 finite RAS `.nii.gz` outputs, and resumes with all stages verified. The synthetic plus workflow composes T1/FLAIR/user-mask/MaSeg/atlas inputs through LesionTOADS, Dura, and CRUISE surfaces/thickness, with restart and tail repair. On the compact Vandy synthetic fixture, a 12-iteration capped Python TOADS run has a final classification hash identical to the selected Java VandyCRUISE result. Committed per-iteration membership arrays localize numeric drift after iteration 3; membership argmax agrees through iteration 10, and the final source classification remains exact even though neither run converges at the cap. Latest full suite on `eb54ed9`: **1,203 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,216.19s**, peak sampled RSS **13.919 GB** under a 24 GB watchdog. Selected TOADS source representation remains an estimate of **3,467–4,867 of 13,651 lines (25–36%)**, with **0 new source lines** this cycle. Whole-workflow Java comparison, intermediate TOADS trajectory/convergence equivalence, isolated-wheel config execution, and real-MRI validation remain open; no human MRI was used. See [12-iteration output and membership evidence](docs/reports/toads-vandy-12-iteration-output-2026-10-06.md), [full-suite result](docs/reports/full-pytest-2026-10-06-vandy-12-membership-trace.md), [compiled base-run evidence](docs/reports/toads-compiled-queue-base-run-2026-10-06.md), [plus workflow evidence](docs/reports/plus-complete-workflow-2026-10-05.md), and the [detailed plan](plan.md).
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
