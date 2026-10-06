# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-06

The restartable brain-resolution synthetic base workflow completes all six stages, converges both TOADS passes (iterations 39 and 25), writes 30 finite RAS `.nii.gz` outputs, and resumes with all stages verified. An isolated wheel runs and resumes all six base stages, verifies 38 manifest output hashes and 30 finite RAS NIfTIs, with each TOADS pass explicitly capped at two iterations and marked nonconverged. The synthetic plus workflow composes T1/FLAIR/user-mask/MaSeg/atlas inputs through LesionTOADS, Dura, and CRUISE surfaces/thickness, with restart and tail repair. The compact Vandy source comparison now preserves Java float32 order for both effective smoothing and rigid-registration energy (`membership * membership * prior * prior`). Python and Java memberships match bit-for-bit through iteration 10; registration transforms match through iteration 11. At the 12-iteration cap, the compact synthetic final labels match voxel-for-voxel and all per-label Dice scores are **1.0**. Neither run is energy-converged, and this remains synthetic fixture evidence rather than whole-workflow or real-MRI validation ([registration energy source-order report](docs/reports/toads-registration-energy-source-order-2026-10-06.md), [Vandy loop report](docs/reports/toads-vandy-source-seeded-loop-diagnostic-2026-10-06.md)). The focused registration, Vandy-loop, and base-restart checks passed **44 tests** on `2619a76`. The alternate sparse TOADS path has a source-trace-matched 343-event thinning queue and a separate 480-event growing queue on folded synthetic inputs; growing includes accepted and rejected topology changes and keeps two membership slots rather than expanding to every class. Full-volume membership/shape validation runs once per thinning pass. Latest full suite on `9c7a800`: **1,209 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,309.68s**, peak sampled RSS **13.363 GB** under a 24 GB watchdog. Selected TOADS source behavior remains estimated at **3,511–4,911 of 13,651 lines (about 26–36%)**; this cycle did not include a full-suite evaluation, so that estimate is unchanged. The approximate route still lacks integrated realistic above-floor behavior, centroid/deviation/membership refresh, registration, full constructor convergence, and output transfer; whole-workflow Java comparison, converged wheel execution, and real-MRI validation remain open. No human MRI was used. See [installed-wheel base report](docs/reports/installed-wheel-base-workflow-capped-2026-10-06.md), [approximate growing evidence](docs/reports/toads-approximate-growing-2026-10-06.md), [latest full-suite result](docs/reports/full-pytest-2026-10-06-approx-growing.md), [approximate thinning evidence](docs/reports/toads-approximate-thinning-pass-2026-10-06.md), [compiled base-run evidence](docs/reports/toads-compiled-queue-base-run-2026-10-06.md), [plus workflow evidence](docs/reports/plus-complete-workflow-2026-10-05.md), and the [detailed plan](plan.md).
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
