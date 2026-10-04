# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-04

The six-stage brain-resolution synthetic base run completes and resumes with all stages verified. Its first TOADS pass reaches the configured 50-iteration cap without its energy threshold; the second converges at iteration 20. The first registered-loop transform is bit-identical to Java at shared inputs, but later growth/thinning and iteration-2 memberships still differ, so whole-TOADS equivalence is not established. This is synthetic integration, not real-MRI validation. LesionTOADS P39zzzzzzbj–bq connect registered initialization through both first-loop branch outcomes and a second synthetic loop transition. The two-loop fixture holds transformed shapes and distance maps fixed; dynamic atlas refresh and convergence remain open. Full suite on `d962c50`: **1,167 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 882.79s**, peak sampled RSS **14.513 GB** under the 24 GB watchdog. Selected Java component evidence plus Python composition does not establish composite Java segmentation or full LesionTOADS equivalence. No human MRI was used. See the [P39bq report](docs/reports/lesion-registered-repeated-iteration-2026-10-04.md), [P39bp report](docs/reports/lesion-registered-first-iteration-continuation-2026-10-04.md), and [plan.md](plan.md).

MaCRUISEplus has an ANTsPy FLAIR-to-T1 registration checkpoint and a `plus-run` CLI that resumes three downstream stages from persisted P39 Dura channels and prealigned maps, producing RAS level sets/thickness and RAS-mm surfaces. Full plus execution, upstream TOADS/LesionTOADS integration, source-workflow comparison, and real-data validation remain open. See [plan.md](plan.md) and the [plus workflow report](docs/reports/plus-file-workflow-2026-10-04.md) for details and reproduction commands.
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
