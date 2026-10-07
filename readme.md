# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-07

The six-stage synthetic base workflow completed all required stages at 172×216×164 and resume verified each stage; that full-run result remains under investigation because its recorded TOADS ratios differ from a new direct source comparison. On the same brain-resolution synthetic N4 input, the current Python TOADS CLI matches selected VandyCRUISE Java through two iterations: both ratios match at float32 precision, the final native-grid labels are voxel-exact (0 changed voxels; Dice 1.0 for all 11 labels), and Python peak sampled RSS was **6.288 GB** under a 16 GB watchdog. This is selected TOADS source evidence, not full-workflow equivalence or 50-iteration convergence. The prior full-run output had 29/30 NIfTI hashes differ from an earlier run and the two first-pass traces are unreconciled; rerunning the complete workflow with implementation fingerprints is the next integration task. The nine-stage raw-input plus workflow completes and naturally converges on the folded 64×80×64 fixture at iteration 6, with finite RAS outputs and verified resume. Compact LesionTOADS growth traces reduced the observed post-fix peak to **0.582 GB** from a preceding **20.221 GB** attempt that exited `-9`. Selected Java comparisons show folded TGDM Dice 1.0, mean surface distance ≤0.003108 mm and Hausdorff ≤0.259171 mm; folded ACE final enhanced-GM and skeleton outputs are voxel-exact. Latest full suite on `7ae1dee`: **1,220 passed, 0 failed, 0 skipped, 11 expected RAS warnings** (1,463.83s, 12.448 GB peak under the 24 GB watchdog); that suite predates the latest ACE/source comparison commits. Next integration work is a clean full base-workflow replay, followed by base/plus source comparisons and real-MRI validation. See the [base rerun report](docs/reports/ace-indexed-heap-brain-resolution-base-run-2026-10-07.md), [full-suite report](docs/reports/full-pytest-2026-10-07-plus-uncrop.md), [folded TGDM comparison](docs/reports/tgdm-folded-cortex-java-comparison-2026-10-07.md), [folded ACE comparison](docs/reports/ace-folded-cortex-java-comparison-2026-10-07.md), and [detailed plan](plan.md).
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
