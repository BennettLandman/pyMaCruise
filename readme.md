# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-03

The source-order-atlas synthetic base workflow completed all six stages on generated `(172,216,164)` inputs. Its manifest verifies RAS native-grid outputs and hashes; a `--resume` run verified and skipped every stage. Peak sampled RSS was 6.427 GB for the run and 0.126 GB for resume, under the 16 GB cap. This is a bounded synthetic integration result, not a converged or source-equivalent result and not real-MRI validation. The [source-atlas report](docs/reports/toads-source-vandy-atlas-comparison-2026-10-03.md) has the runnable command, output hashes, and stage evidence.

Both TOADS passes reached their two-iteration caps (terminal ratios 1.10402799 and 0.571384311). The source-multigrid port now matches the selected Java implementation on the full `(163,207,156)` crop: maximum component error is `6.77e-6`, below the `5e-5` gate. The base workflow recomputed core and refinement from this solver through its resume path; all stage checks passed. TOADS remains capped at two iterations per pass, so the base result is still nonconverged and is not whole-workflow Java equivalence. The compact Java/Python TOADS comparison matches startup and the first cycle closely; after 12 iterations, final label Dice is at least 0.9827 on all 11 labels, although energy ratios diverge after iteration 2. These selected synthetic comparisons do not establish whole-workflow numerical acceptance.

The source-order loader change passed its focused regression command (**14 passed in 9.64s**, code commit `4b57532`). After the multigrid port, focused GVF/plus-GVF/CLI/core/restart tests passed (**19 passed in 77.68s**); the full suite has not yet been rerun on that change. The latest full-suite evidence is **1,138 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 793.64s** on `ad34321`, before the loader change; the full suite has not been rerun. No real-MRI validation is recorded. See the [full-suite report](docs/reports/full-pytest-2026-10-03-brain-resolution.md).

**Next milestone:** diagnose source-order TOADS convergence at brain resolution, then advance the plus workflow with validated base-stage settings. Keep synthetic integration, Java source comparisons, whole-workflow numerical acceptance, and real-data validation distinct. See [plan.md](plan.md) for row-level progress.

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
