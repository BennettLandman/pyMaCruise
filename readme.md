# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-11

The restartable base workflow completed all eight stages on the supplied human T1, brain mask, and MaSeg after regridding oblique inputs to a shared orthogonal RAS grid. The run used the selected V3 exposed settings (including five Dura iterations and TGDM initial/inner/central/outer counts 6/4/7/7); a same-settings resume verified all eight stages. It follows the chosen native-RAS path with a caller-supplied mask rather than historical BET/SPECTRE preprocessing. Raw MaSeg support remains below target: 82.34% GM and 88.48% WM (containment, not Dice), and refinement leaves supported counts unchanged. The three watertight meshes have 143 components each; direct Java surface extraction matches Python on the fixed fields within the distance limits. Some lower-level Dura parameter mapping remains under review, and no historical full-pipeline output is available, so full workflow equivalence and clinical validity remain unestablished. Human inputs/outputs remain under ignored `build/`. See the [full-chain human report](docs/reports/full-chain-human-base-2026-10-11.md), [restartable-tail report](docs/reports/restartable-human-base-tail-2026-10-11.md), [Java surface comparison](docs/reports/in-vivo-full-human-jist-surface-extraction-2026-10-11.md), and [detailed plan](plan.md).

The latest full repository suite was run on commit `3bb0634`: **1,512 passed, 1 failed, 0 skipped, 11 warnings** in 1,137.75 s, with 12.730 GB peak sampled RSS. Its failure was a stale plus-resume test expectation; a focused correction passed, but the full suite has not been rerun since. Current focused workflow validation: **12 passed in 29.48 s**; no full suite was run for this cycle. The source-method estimate remains **3,511–4,911 of 13,651 selected TOADS lines (26–36%)**; this does not measure workflow completion. Current blockers are inadequate human ribbon support and unavailable full-pipeline reference outputs, broader standalone TOADS/LesionTOADS coverage, plus-workflow quality/parity, CLI completion, and paper annotations. Synthetic evidence remains distinct from Java component comparisons and real-MRI integration. See [plan.md](plan.md), [folded-cortex validation](docs/reports/labeled-folded-cortex-validation-2026-10-09.md), and [full-suite evidence](docs/reports/full-pytest-2026-10-10-dura03.md).

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
