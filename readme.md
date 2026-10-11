# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-11

The full human BASE chain is being rebuilt from the hash-verified 1-mm RAS T1, user brain mask, and MaSeg with current source code, including both natural-convergence TOADS passes, cleaned MaSeg CSF boundary, and the previously selected parameters. The run is under a 16 GB RSS watchdog; it has reached the first TOADS pass. Earlier boundary-enabled tail metrics—82.3790% GM and 88.4888% WM containment, not Dice, with 143 components per surface—are historical until this full rebuild completes. Held-fixed Java/Python replay after the ACE FMM correction matches enhanced GM, cleaned skeleton, and thinned skeleton across all 3,195,483 voxels; this does not establish whole-workflow equivalence. See the [corrected human-tail report](docs/reports/in-vivo-base-tail-with-v3-boundary-2026-10-11.md), [ACE correction report](docs/reports/in-vivo-human-ace-ordered-fmm-correction-2026-10-11.md), [ACE stage diagnosis](docs/reports/in-vivo-human-ace-stage-localization-2026-10-11.md), [boundary integration evidence](docs/reports/base-ace-boundary-integration-2026-10-11.md), and [detailed plan](plan.md).

The TOADS constructor-border and loop-carried skeleton state follows selected Java behavior. Full suite on implementation commit `9cade54`: **1,531 passed, 0 failed, 0 skipped, 11 expected RAS warnings** in 801.22 s; peak sampled RSS was 13.450 GB under the 16 GB watchdog. The selected TOADS source-line estimate remains **3,511–4,911 of 13,651 lines (26–36%)**; no full suite has been run for the newer P41z changes. P41z iteration-6 thinning matches Java across all **267,638 queue events**, its 13-class statistics refresh is bitwise exact, and the full held-fixed membership update now matches every float32 output bit after correcting a raw-lesion-energy omission in the Java harness. The plus run remains nonterminal at energy ratio 0.131597 vs 0.0001 and fails ribbon-quality checks. Synthetic component matches and held-fixed human source comparisons remain distinct from complete workflow validation. See the [latest full-suite report](docs/reports/full-pytest-2026-10-11-toads-skeleton-carry.md), [P41z membership checkpoint](docs/reports/p41z-iteration6-membership-carry-2026-10-11.md), [P41z thinning](docs/reports/p41z-iteration6-thinning-java-2026-10-11.md), [P41z statistics](docs/reports/p41z-iteration6-statistics-java-2026-10-11.md), and [detailed plan](plan.md).

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
