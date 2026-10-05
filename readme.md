# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-05

The six-stage brain-resolution synthetic base run completes and resumes with all stages verified. Its first TOADS pass reaches the configured 50-iteration cap without its energy threshold; the second converges at iteration 20. The first registered-loop transform is bit-identical to Java at shared inputs, but later growth/thinning and iteration-2 memberships still differ, so whole-TOADS equivalence is not established. LesionTOADS P39zzzzzzbj–bv connect registered initialization through synthetic loop updates, refreshed atlas shapes and distance maps, and stop exits; P41o/p add restartable registered continuation and terminal classification/subtype outputs. P41q now adds `pymacruise plus-full-run`, from raw T1/FLAIR/user mask and atlas descriptor through preprocessing to a terminal or continuing LesionTOADS checkpoint, with stage hashes, resume, and corrupt-tail repair. This is synthetic integration, not a Java composite match or real-MRI validation; Dura-to-CRUISE surface processing remains downstream. Latest full suite on tested commit `7dc8a67` passed **1,193 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 995.58s**, with peak sampled RSS **14.609 GB** under the 24 GB watchdog. The selected base-TOADS estimate remains 3,464–4,864 of 13,651 lines (25–36%); this workflow-only slice adds 0 mapped Java method lines, and LesionTOADS lines remain outside that denominator. No human MRI was used. See the [raw-input plus workflow report](docs/reports/plus-full-workflow-2026-10-05.md), [terminal checkpoint report](docs/reports/plus-terminal-lesion-checkpoint-2026-10-05.md), [base-run report](docs/reports/base-workflow-brain-resolution-converged-2026-10-04.md), [latest full-suite report](docs/reports/full-pytest-2026-10-05-plus-full-workflow.md), and [plan.md](plan.md).

The raw-input plus workflow now reaches terminal LesionTOADS outputs, with a separate `plus-run` continuation available from persisted P39 Dura channels and prealigned maps. That downstream command is not yet wired to the terminal outputs from `plus-full-run`; completing that handoff and Dura-to-CRUISE surface workflow, then source comparison and real-data validation, remains open. See the [plus workflow report](docs/reports/plus-full-workflow-2026-10-05.md), [plan.md](plan.md), and the [CLI input/output reference](docs/cli-help.md).
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
