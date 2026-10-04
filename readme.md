# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-04

The six-stage brain-resolution synthetic base run completes and resumes with all stages verified. Its first TOADS pass reaches the configured 50-iteration cap without its energy threshold; the second converges at iteration 20. A selected-Java/Python synthetic trace now localizes the earliest observed TOADS split after iteration-1 memberships: six registration parameters differ by about 0.0002, followed by 55/162,180 post-growth labels and 843/648,720 closest-label entries. The later thinning boundary differs at 42 voxels (minimum Dice 0.9956), and iteration-2 memberships diverge locally (maximum absolute difference 0.394). The cause remains open. The base run is synthetic integration, not whole-workflow source equivalence or real-MRI validation. LesionTOADS has source-pinned component and fixed-state membership evidence, but not a complete iterative workflow.

MaCRUISEplus now has an ANTsPy FLAIR-to-T1 registration checkpoint and a `plus-run` CLI that resumes three downstream stages from persisted P39 Dura channels and prealigned maps, producing RAS level sets/thickness and RAS-mm surfaces. The repository suite passed **1,158 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 844.58s**, peak RSS **14.536 GB** under the 24 GB watchdog on implementation commit `4ae3b55`. Later narrow changes to interrupted-run handling and T1-grid validation have focused test evidence on `371a807` and `8307df3`; the actual CLI resume test passed on `e741c65`. Those later changes have not had another full-suite run. Full plus execution, upstream TOADS/LesionTOADS integration, source-workflow comparison, and real-data validation remain open. See [plan.md](plan.md) and the [plus workflow report](docs/reports/plus-file-workflow-2026-10-04.md) for details and reproduction commands.
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
