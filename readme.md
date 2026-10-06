# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-06

The restartable brain-resolution synthetic base workflow completes all six stages, converges both TOADS passes (iterations 39 and 25), writes 30 finite RAS `.nii.gz` outputs, and resumes with all stages verified. An isolated wheel now also runs and resumes all six base stages, verifies 38 manifest output hashes and 30 finite RAS NIfTIs, with each TOADS pass explicitly capped at two iterations and marked nonconverged. The synthetic plus workflow composes T1/FLAIR/user-mask/MaSeg/atlas inputs through LesionTOADS, Dura, and CRUISE surfaces/thickness, with restart and tail repair. On the compact Vandy synthetic fixture, a 12-iteration capped Python TOADS run has a final classification hash identical to the selected Java VandyCRUISE result, but neither run converges at the cap; small membership differences start at the first exact update and amplify later. Java-state-seeded first and fourth updates match within one float32 ULP, and the initial source/Python memberships, startup standard deviations, closest-label map, shape samples, and first thinning labels match exactly ([diagnostic](docs/reports/toads-vandy-source-seeded-loop-diagnostic-2026-10-06.md)). The alternate sparse TOADS path now has a source-trace-matched 343-event thinning queue and a separate 480-event growing queue on folded synthetic inputs; growing includes accepted and rejected topology changes and keeps two membership slots rather than expanding to every class. Full-volume membership/shape validation runs once per thinning pass. Latest full suite on `9c7a800`: **1,209 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,309.68s**, peak sampled RSS **13.363 GB** under a 24 GB watchdog. Selected TOADS source behavior is estimated at **3,511–4,911 of 13,651 lines (about 26–36%)**; this cycle adds about **7 selected Java lines** for approximate growing, after the preceding thinning queue’s 37. The approximate route still lacks integrated realistic above-floor behavior, centroid/deviation/membership refresh, registration, full constructor convergence, and output transfer; whole-workflow Java comparison, intermediate TOADS trajectory equivalence, converged wheel execution, and real-MRI validation remain open. No human MRI was used. See [installed-wheel base report](docs/reports/installed-wheel-base-workflow-capped-2026-10-06.md), [approximate growing evidence](docs/reports/toads-approximate-growing-2026-10-06.md), [latest full-suite result](docs/reports/full-pytest-2026-10-06-approx-growing.md), [approximate thinning pass evidence](docs/reports/toads-approximate-thinning-pass-2026-10-06.md), [full-suite result](docs/reports/full-pytest-2026-10-06-approx-thinning-validation-boundary.md), [compiled base-run evidence](docs/reports/toads-compiled-queue-base-run-2026-10-06.md), [plus workflow evidence](docs/reports/plus-complete-workflow-2026-10-05.md), and the [detailed plan](plan.md).
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
