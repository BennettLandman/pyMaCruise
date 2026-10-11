# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-11

A source audit found that the earlier eight-stage human base run selected ACE mode 2 but omitted its cleaned MaSeg CSF boundary, leaving that ACE behavior inactive. The previous human support and mesh results are provisional diagnostics and are not accepted as a complete V3 base run. The boundary path is now connected in both restartable base workflows; **18 synthetic integration tests passed in 15.86 s** on commit `0d416ff`. A corrected human tail rerun from the verified second-pass TOADS checkpoint is underway. It preserves the user-selected native-RAS/caller-mask path and does not run historical BET/SPECTRE preprocessing. See the [boundary integration report](docs/reports/base-ace-boundary-integration-2026-10-11.md) and [detailed plan](plan.md).

The latest full repository suite remains the historical run on `3bb0634`: **1,512 passed, 1 failed, 0 skipped, 11 warnings** in 1,137.75 s, with 12.730 GB peak sampled RSS; a focused correction addressed its stale plus-resume test expectation, but the full suite has not been rerun since. No full suite was run for this boundary integration, so the selected TOADS source-method estimate remains **3,511–4,911 of 13,651 (26–36%)**, with zero mapped Java lines added. Remaining work includes human ribbon-quality acceptance, Dura parameter interpretation, full-pipeline reference comparison, plus-workflow quality/parity, broader TOADS/LesionTOADS behavior, standalone CLIs, and paper annotations. Synthetic, fixed-component source comparison, and human integration evidence remain separate categories.

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
