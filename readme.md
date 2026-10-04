# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-04

The restartable brain-resolution synthetic base workflow has run through all six stages, and a same-parameter resume verified and skipped every stage. Outputs are finite RAS NIfTI on the native grid; peak sampled RSS was 8.308 GB under a 16 GB watchdog. The first TOADS pass hit its 50-iteration cap without reaching the energy threshold, while the second converged at iteration 20. The source V3 wrapper itself configures a `1e-4` threshold and 50-iteration cap, so this stop is allowed by the selected source settings; it is not by itself evidence of a Python defect. The manifest records `status=complete`, `converged=false`, and `source_equivalence_claim=false`. The brain-resolution inputs are committed under `tests/data/synthetic/base_workflow_brain_resolution/`. Full-route Java output comparison and clinical validation remain open. Details and reproduction commands are in the [brain-resolution base-run report](docs/reports/base-workflow-brain-resolution-converged-2026-10-04.md).

The exact polynomial field refit accepts descriptor-typed 13-class LesionTOADS data and is connected to a synthetic source-matched adaptive two-channel transition. Active raw energy dispatch for the bundled 13-class descriptor matches the pinned Java harness bit-for-bit for OBJECT, LESION, MASK, and optimized-CSF fallback terms. The descriptor’s fixed-term membership update also matches all class, lesion, and denominator bits; the combined raw-energy-to-membership transition matches the pinned one-voxel Java case. A statistics-to-membership adapter exercises optional FLAIR lesion weighting on synthetic fixed state. Its focused suite passes **8 tests** and the full suite passes **1,153 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 828.21s** on tested commit `936ba20`; peak RSS was not sampled. P39zzzzzzbg is Done for this fixed-state composition; statistics-refreshed iterative LesionTOADS execution remains open. Plus P41a connects generated P39 Dura channels through transfer/cropping, topology, ACE, GVF, TGDM, uncropping, surfaces, and thickness on committed synthetic atlas/TOADS inputs. Its full suite passes **1,154 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 834.21s** on tested commit `6779753`, with peak sampled RSS of 15.381 GB under a 24 GB watchdog. Atlas/TOADS inputs are precomputed; complete LesionTOADS and plus restart/manifest remain open. P41b adds a restartable ANTsPy FLAIR-to-T1 checkpoint and `pymacruise plus-register-flair` CLI, with focused synthetic validation; NiftyReg/Aladin equivalence is unvalidated. These are synthetic component/integration results, not workflow equivalence or real-MRI validation. Remaining base-workflow evidence includes source-matched Java/Python comparison of TOADS loop outputs and stage-wise workflow outputs. The plus workflow, real-MRI validation, and published-paper documentation/annotations remain open. See [plan.md](plan.md) for the detailed steps and current blockers.

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
