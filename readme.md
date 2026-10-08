# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-07

The corrected six-stage synthetic base workflow completes and resumes on the `(172,216,164)` brain-resolution fixture. All 30 manifest NIfTIs are finite and RAS; both TOADS passes converge (12 and 34 iterations), all six stage checks pass, and resume skips all six stages. Selected Java comparisons now match the first pass through its natural stop, the corrected Dura stage within agreed tolerances (22 mask-boundary differences; Dice 0.9999900; symmetric mean distance 0.000205 mm; Hausdorff 0.9375 mm), and the second pass through natural convergence (34 float32 energy ratios and all final labels exact). Peak sampled RSS was **8.106 GB** under the 16 GB watchdog; resume used **0.127 GB**. Focused checks on implementation commit `9eb84bc`: TGDM **24 passed in 14.85s**, restartable base workflow **3 passed in 41.13s**. The full suite on tested commit `c17d2c4` passes **1,230 tests, 0 failures, 0 skips, 11 expected RAS warnings** in **1,607.32s**, with **13.364 GB** peak RSS under a 16 GB watchdog; see [full-suite report](docs/reports/full-pytest-2026-10-07-plus-resume-versioning.md). Evidence and exact reproduction details: [corrected base workflow report](docs/reports/base-workflow-after-tgdm-topology-fix-2026-10-07.md). These are synthetic stage and selected-source comparisons, not whole-workflow source equivalence or real-MRI validation. Base refinement/source comparison and plus source parity remain open. The nine-stage raw-input plus workflow was freshly rerun after the TGDM correction on the 64×80×64 folded fixture: all stages completed, LesionTOADS converged at iteration 6, and resume verified all nine stages. Peak RSS was **0.810 GB** under a 16 GB watchdog; all 16 scientific final products match the prior same-input run. This is Python integration evidence, not plus source equivalence; see [plus rerun report](docs/reports/plus-after-tgdm-correction-2026-10-07.md). Implementation `caf876d` also fingerprints plus code for restart: changed code runs in a new preserved stage generation. Two focused restart regressions pass (**1 in 102.89s**, **1 in 79.55s**; maximum sampled RSS **0.582 GB**); these are included in the full suite on `c17d2c4`. See the [brain-resolution Dura comparison](docs/reports/dura-brain-resolution-java-comparison-2026-10-07.md), [CRUISE core comparison](docs/reports/tgdm-core-brain-resolution-java-comparison-2026-10-07.md), [outer TGDM source trace](docs/reports/tgdm-outer-shared-central-java-comparison-2026-10-07.md), [Dura wrapper comparison](docs/reports/dura-wrapper-41-java-comparison-2026-10-07.md), [MaACE brain-resolution comparison](docs/reports/ace-brain-resolution-java-comparison-2026-10-07.md), [latest full-suite report](docs/reports/full-pytest-2026-10-07-plus-resume-versioning.md), and [detailed plan](plan.md).
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
