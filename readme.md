# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-08

The full-registration brain-resolution synthetic base run completed all seven stages, with both TOADS passes naturally converging (11 and 48 iterations). An independent audit verified 46/46 output hashes and 31 finite RAS NIfTIs; the 16 GB watchdog recorded an 8.514 GB peak. Same-settings resume verified and skipped all seven stages under a 4 GB watchdog. This end-to-end integration run took 3:08:33; see the [full-registration base workflow report](docs/reports/base-workflow-natural-registration-2026-10-08.md).

The folded synthetic plus workflow completes and resumes all nine stages; LesionTOADS converges at iteration 6. Its latest audit verified 75 manifest-referenced hashes and 62 finite RAS NIfTIs. Selected-source checks include direct VandyCRUISE rigid-pyramid comparisons and the generated JIST connectivity-18 case table. The latest full suite on implementation tree `5a0adcf` passed **1,502 tests, 0 failed, 0 skipped, and 11 expected RAS warnings** in 1,873.83 s at 12.442 GB peak RSS; code did not change after that tested tree. The selected TOADS source-method estimate remains 3,511–4,911 of 13,651 lines (about 26–36%); this workflow cycle mapped 0 additional source lines.

These results establish synthetic integration and selected-source behavior, not whole-workflow MATLAB/JIST equivalence or real-MRI validation. Remaining work includes source-composite comparisons, MaSeg refinement provenance and source comparison, plus-workflow parity, paper-based annotations, and real-data validation. The available human T1 still needs its matching user-provided brain mask. See the [base report](docs/reports/base-workflow-natural-registration-2026-10-08.md), [plus report](docs/reports/plus-workflow-current-pyramid-source-2026-10-08.md), and [full-suite report](docs/reports/full-pytest-2026-10-08-bundled-atlas-startup.md).

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
