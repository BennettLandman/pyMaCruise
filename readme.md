# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-08

The six-stage brain-resolution synthetic base run was rebuilt after the latest thickness changes and completed successfully: six stage checks passed, both TOADS passes converged at the previously matched Java stop points (12 and 34 iterations), all manifest NIfTIs are finite RAS, and 16 final products were recorded. Peak sampled RSS was **8.190 GB** under a 16 GB watchdog. A same-settings resume verified and skipped all six stages under a 4 GB watchdog at **0.123 GB**. Focused restart tests: **3 passed in 41.40s** on `b500127`; no full suite was run in this cycle. Reproduction and evidence: [2026-10-08 base restart report](docs/reports/base-workflow-restart-after-thickness-change-2026-10-08.md).

This validates synthetic integration and restart behavior, not whole-workflow source equivalence or real-MRI performance. The base workflow emits tissue/ACE/GVF images, level sets, and refined labels; it does not generate meshes or thickness maps. The selected JIST thickness wrapper separately matches its Java reference on held-fixed meshes. Open work includes whole-workflow and refinement comparisons, JIST-versus-Lewiner topology differences, plus source parity, remaining selected TOADS/LesionTOADS coverage, paper annotations, and real-MRI validation. The provided Singularity output directory is empty and the available human T1 lacks a matching brain mask. The last full suite remains **1,234 passed** on `e774747`; see [full-suite report](docs/reports/full-pytest-2026-10-07-lesion-objective-source-match.md).
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
