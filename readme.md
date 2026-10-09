# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-09

Synthetic base and plus workflows are restartable and produce their documented RAS NIfTI outputs; the current brain-resolution base run stayed under the 16 GB target and the plus run under 0.8 GB. Selected-source Java comparisons now cover TOADS/Dura handoffs, MaACE, TGDM, JIST surfaces, and thickness at documented held-fixed inputs. The labeled 65³ folded-cortex case uses the full source-selected TGDM iteration profile and explicit MaSeg-like WM/GM/CSF labels. Its production JIST meshes pass topology, nesting, label-support, and ribbon-quality checks; triangle closest-point thickness averages 4.673 mm and level-set subtraction averages 4.643 mm. All Java TGDM final-surface gates pass on this case. These are synthetic/stage-level results, not full composite equivalence or human MRI validation.

Latest full suite: **1,505 passed, 0 failed, 0 skipped, 11 expected RAS warnings** in 1,831.87 s, with 12.464 GB peak RSS under a 16 GB watchdog, on tested implementation commit `eb04598`. The subsequent validation-runner/report commit `7111e4a` has direct bounded-run and Java-comparison evidence but is not part of that full-suite run. The coarse selected TOADS source-method estimate remains **3,511–4,911 of 13,651 lines (26–36%)**; it is not a workflow-completion metric. Remaining blockers include whole-workflow source parity, broader TOADS/LesionTOADS behavior, standalone CLIs, paper annotations, and real-MRI validation (the provided T1 still needs its matching brain mask). See [plan.md](plan.md), [folded-cortex validation](docs/reports/labeled-folded-cortex-validation-2026-10-09.md), and [full-suite evidence](docs/reports/full-pytest-2026-10-09-labeled-folded-cortex.md).

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

The runtime stack is NiBabel, NumPy, and SciPy; pytest, build, and scikit-image are development dependencies, and `antspyx` is optional for ANTs-backed operations. The license/attribution inventory and follow-up TTO work are in [license_info.md](license_info.md). All atlas images needed at runtime are bundled in this repository as regular Git files.

## Inputs and outputs

The command-line input, assumption, and output contract is documented in [docs/cli-help.md](docs/cli-help.md). Coordinate and NIfTI rules are in [docs/image-coordinate-conventions.md](docs/image-coordinate-conventions.md). Numerical acceptance thresholds are in [docs/validation-tolerances.md](docs/validation-tolerances.md).
