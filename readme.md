# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-05

`pymacruise plus-full-run` composes raw T1/FLAIR/user-mask/MaSeg/atlas inputs through a restartable synthetic run ending in CRUISE level sets, surfaces, and thickness outputs. The folded-cortex case resumes and repairs corrupted tail outputs; surface thickness is finite, while volume-thickness reports unavailable when there are no crossing voxels. `pymacruise workflow-run --config workflow.json [--resume]` runs and resumes synthetic base and plus workflows from the source checkout. LesionTOADS registered-loop correction fields now persist as hashed RAS `.nii.gz`, reject uncompressed field inputs, and carry adaptive degree through resumed loop segments; both the single-segment CLI and restartable loop CLI expose the feature. Latest full suite on implementation commit `c8a3303`: **1,202 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,187.29s**, peak sampled RSS **13.539 GB** under the 24 GB watchdog. The selected TOADS source estimate is 3,467–4,867 of 13,651 lines (25–36%). The first synthetic growth, closest-label state, and next thinning now match the selected Java reference exactly after carrying thinning priorities into growth; iteration-2 memberships differ by at most `4.31e-6`. The brain-resolution base run's first TOADS pass reaches its 50-iteration cap without meeting its energy threshold; the second converges at iteration 20, and later TOADS outputs still differ from Java. The isolated wheel exposes the command and resolves bundled assets, but full workflow execution from the wheel is open. No human MRI was used. See the [TOADS ordering-carry report](docs/reports/toads-ordering-carry-2026-10-05.md), the [P39 field-checkpoint report](docs/reports/lesion-field-checkpoints-2026-10-05.md), [latest full-suite report](docs/reports/full-pytest-2026-10-05-toads-ordering-carry.md), [P41t synthetic workflow report](docs/reports/plus-complete-workflow-2026-10-05.md), [plan.md](plan.md), and [CLI input/output reference](docs/cli-help.md).
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
