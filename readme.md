# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-02

**The restartable, single-command synthetic base MaCRUISE workflow is still being assembled.** Existing integration tests compose N4 → first TOADS → Dura → second TOADS → optional MaSeg cleanup → ACE → GVF → nested TGDM and produce native-grid RAS level sets, but do not yet provide verified stage restart, a complete run manifest, or convergence evidence. Compact TOADS runs hit their iteration cap; the Java and Python stopping rules match, so the threshold is not being relaxed.

- **Runnable diagnostic:** `./.venv/bin/python scripts/trace_synthetic_toads_loop.py --output build/benchmarks/synthetic-toads-10.nii.gz --maximum-iterations 10 --registration-levels 1 --registration-iterations 1 --loop-registration-iterations 1 --padding 2 --relation-limit 1`. It writes uint8 RAS labels 1–11; the 10-cycle run ended at ratio `0.545753181`, iteration cap reached, with output SHA-256 `779c64ba73059d716fd65fcbbec3193a2e0bd5701cd687e223c530650979d764`. This is a repeatable synthetic runtime/stop trace, not Java output equivalence or the full base run. See [the characterization report](docs/reports/synthetic-toads-loop-characterization-2026-10-02.md).
- **Latest full suite:** `./.venv/bin/pytest -q` passed **1,110 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 1,305.91 seconds** on tested commit `8b4a13a9600e94b7a33c60cbbf6f57eafda9b876`; see [the full-suite report](docs/reports/full-pytest-2026-10-02-folded-integration.md). Later approximate-TOADS and LesionTOADS commits have focused evidence only and were not part of that suite.
- **Newly covered source behavior:** Java-shaped nonfinite TOADS registration-trial rejection; TOADS folded-cortex thinning/update/growing/relation integration; LesionTOADS atlas startup/first thinning, output CLI, and constructor state. Separate source comparisons and Python-only integration evidence are identified in their reports.
- **Remaining blockers:** complete and restartable base orchestration with required images, stage checks, and manifest; characterize TOADS behavior without misreporting cap termination as convergence; then establish source/output comparisons and advance MaCRUISEplus. Human MRI validation remains deferred. Historical SPECTRE transfer is intentionally omitted under native-T1-RAS operation with a user-supplied T1-grid brain mask.

Checklist percentages and aggregate test counts are bookkeeping, not measures of workflow completion. See [plan.md](plan.md) for step evidence and the next integration work.
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
