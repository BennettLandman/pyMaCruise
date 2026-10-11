# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-10

Synthetic base and plus workflows are restartable and emit RAS NIfTI products. The committed 65³ labeled folded-cortex TGDM case passes explicit ribbon-quality checks and all Java surface gates (Dice 1.0; maximum Hausdorff 0.070323 mm); a held-fixed 73×91×63 TGDM stage matches Java across five surfaces (zero sign differences, Dice 1.0, maximum Hausdorff 0.0000221 mm). The supplied human MRI completed the corrected source-profile base chain, including second TOADS convergence at 22 iterations and surface/thickness products. In a direct Java replay, pre-evolution TOADS classifications match exactly; after two iterations, 7,057 voxels differ but every label Dice is ≥0.980. At natural convergence, Java stops at iterations 23/24 while Python stops at 22, and labels 6/16 fall below the 0.97 Dice threshold. GM/WM coverage remains below the 97% ribbon QA check, and each mesh has 209 components. No historical MaCRUISE outputs were available for final surface comparison. This is real-data integration and partial component-level source evidence, not whole-pipeline equivalence or clinical validity. Human images/products remain under ignored `build/invivo-base-2026-10-10-ras-grid/`; see the [in-vivo report](docs/reports/in-vivo-base-workflow-smoke-2026-10-10.md), [TOADS replay report](docs/reports/in-vivo-toads-java-localization-2026-10-10.md), and [plan](plan.md).

Latest full suite: **1,512 passed, 1 failed, 0 skipped, 11 warnings** in 1,137.75 s, with 12.730 GB peak sampled RSS under a 24 GB watchdog, on implementation commit `3bb0634`. The failure was a stale plus-resume test expectation; the corrected focused test passes (1 passed in 41.37 s), but the full suite has not been rerun after that test-only correction. The selected TOADS source-method estimate remains **3,511–4,911 of 13,651 lines (26–36%)**, with zero newly mapped source lines in this cycle; this is not a workflow-completion metric. Remaining blockers include whole-workflow source parity, TOADS natural-convergence differences, broader TOADS/LesionTOADS behavior, standalone CLIs, paper annotations, and in-vivo ribbon-quality/reference comparisons. See [plan.md](plan.md), [folded-cortex validation](docs/reports/labeled-folded-cortex-validation-2026-10-09.md), [TOADS replay report](docs/reports/in-vivo-toads-java-localization-2026-10-10.md), and [suite and focused-test evidence](docs/reports/full-pytest-2026-10-10-dura03.md).

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
