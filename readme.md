# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-04

The six-stage brain-resolution synthetic base run completes and resumes with all stages verified. Its first TOADS pass reaches the configured 50-iteration cap without its energy threshold; the second converges at iteration 20. A selected-Java/Python synthetic trace localizes the earliest observed TOADS split after iteration-1 memberships: six registration parameters differ by about 0.0002, followed by 55/162,180 post-growth labels and 843/648,720 closest-label entries. The later thinning boundary differs at 42 voxels (minimum Dice 0.9956), and iteration-2 memberships diverge locally (maximum absolute difference 0.394). Iteration-2 centroids remain within 9.54e-7; five class standard deviations differ on the separate label states, but the Python variance kernel matches Java within 2e-7 when given identical Java state. The registration split is traced to a 1.93e-8 post-refresh transform difference: with the same Java transform, Python reproduces the Java priors and baseline energy exactly; the tiny transform change shifts energy by 1.79e-7 and flips the strict line-search branch. The upstream rounding source remains open. The base run is synthetic integration, not whole-workflow source equivalence or real-MRI validation. LesionTOADS now composes a 13-class transition through next thinning, statistics refresh, GM membership, rigid atlas registration, and the wrapper's stop/growth/classification branch; the assembled transform and outputs match the manually composed Java-pinned synthetic path. Spatial/distance refresh and correction-field fitting remain external to that transition, and the full iterative workflow remains open.

MaCRUISEplus has an ANTsPy FLAIR-to-T1 registration checkpoint and a `plus-run` CLI that resumes three downstream stages from persisted P39 Dura channels and prealigned maps, producing RAS level sets/thickness and RAS-mm surfaces. The repository suite passed **1,159 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 862.69s**, peak sampled RSS **14.584 GB** under the 24 GB watchdog on implementation commit `7305647`; see the [full-suite report](docs/reports/full-pytest-2026-10-04-lesion-registered-iteration.md). Full plus execution, upstream TOADS/LesionTOADS integration, source-workflow comparison, and real-data validation remain open. See [plan.md](plan.md) and the [plus workflow report](docs/reports/plus-file-workflow-2026-10-04.md) for details and reproduction commands.
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
