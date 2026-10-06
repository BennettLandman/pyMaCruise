# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-05

`pymacruise base-run` now completes a restartable six-stage brain-resolution synthetic run with both TOADS passes converged (iterations 39 and 25), all stage checks passed, finite RAS native-grid final outputs, and all stages verified/skipped on resume. The run peaked at 6.887 GB RSS under the 16 GB watchdog. This is synthetic integration, not whole-workflow Java/MATLAB equivalence or human-MRI validation. `pymacruise plus-full-run` composes raw T1/FLAIR/user-mask/MaSeg/atlas inputs through a restartable synthetic run ending in CRUISE level sets, surfaces, and thickness outputs; the folded-cortex case resumes and repairs corrupted tail outputs. The source checkout's `workflow-run --config workflow.json [--resume]` runs and resumes synthetic base and plus workflows. The installed wheel's `plus-full-run` completed and resumed all nine synthetic plus stages outside the source checkout; config-driven `workflow-run` from a wheel remains open. LesionTOADS registered-loop correction fields persist as hashed RAS `.nii.gz` and carry adaptive degree through resumed loop segments; both single-segment and restartable loop CLIs expose the feature. Latest full suite on implementation commit `c8a3303`: **1,202 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,187.29s**, peak sampled RSS **13.539 GB** under the 24 GB watchdog. The selected TOADS source estimate is 3,467–4,867 of 13,651 lines (25–36%). The selected growth, closest-label state, and next thinning match the Java reference exactly; iteration-2 membership maximum difference is `4.31e-6`. Whole-workflow source comparison and real-MRI validation remain open. No human MRI was used. See the [brain-resolution base-run report](docs/reports/base-workflow-ordering-carry-2026-10-05.md), [installed-wheel report](docs/reports/installed-wheel-plus-workflow-2026-10-05.md), [TOADS ordering-carry report](docs/reports/toads-ordering-carry-2026-10-05.md), [P39 field-checkpoint report](docs/reports/lesion-field-checkpoints-2026-10-05.md), [latest full-suite report](docs/reports/full-pytest-2026-10-05-toads-ordering-carry.md), [P41t synthetic workflow report](docs/reports/plus-complete-workflow-2026-10-05.md), [plan.md](plan.md), and [CLI input/output reference](docs/cli-help.md).
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
