# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-05

The six-stage brain-resolution synthetic base run completes and resumes with all stages verified. Its first TOADS pass reaches the configured 50-iteration cap without its energy threshold; the second converges at iteration 20. The first registered-loop transform is bit-identical to Java at shared inputs, but later growth/thinning and iteration-2 memberships still differ, so whole-TOADS equivalence is not established. The restartable `pymacruise plus-full-run` now covers raw T1/FLAIR/user-mask inputs through terminal or continuing LesionTOADS output, four selected Dura-removal channels on the registered recrop grid, and GM/WM-fill/WM-mask transferred to native T1 RAS geometry. MaSeg-derived masks and downstream CRUISE surface/thickness stages remain open. This is synthetic integration, not Java composite equivalence or real-MRI validation. Latest full suite on tested commit `2c32fd8` passed **1,193 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 1,026.32s**, with peak sampled RSS **14.479 GB** under the 24 GB watchdog. The selected base-TOADS estimate remains 3,464–4,864 of 13,651 lines (25–36%); this workflow composition adds 0 mapped Java method lines, and LesionTOADS lines remain outside that denominator. No human MRI was used. See the [native transfer report](docs/reports/native-dura-transfer-2026-10-05.md), [Dura handoff report](docs/reports/lesion-dura-output-checkpoint-2026-10-05.md), [raw-input plus workflow report](docs/reports/plus-full-workflow-2026-10-05.md), [latest full-suite report](docs/reports/full-pytest-2026-10-05-native-dura-transfer.md), and [plan.md](plan.md).

`plus-full-run` now creates selected Dura-removal channels after its terminal LesionTOADS checkpoint and transfers the three CRUISE inputs to native T1 RAS geometry. MaSeg-derived tissue masks and composition through the downstream surface/thickness outputs remain open. The physical-affine transfer adapter has not been compared against JIST for arbitrary differing grids. See the [native transfer report](docs/reports/native-dura-transfer-2026-10-05.md), [plan.md](plan.md), and the [CLI input/output reference](docs/cli-help.md).
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
