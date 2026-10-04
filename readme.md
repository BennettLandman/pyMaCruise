# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-04

The source-order-atlas synthetic base workflow completed all six stages on generated `(172,216,164)` inputs. Its manifest verifies RAS native-grid outputs and hashes; a `--resume` run verified and skipped every stage. Peak sampled RSS was 6.427 GB for the run and 0.126 GB for resume, under the 16 GB cap. That base route used deliberately capped two-iteration TOADS passes. Separately, the standalone source-order brain-resolution TOADS CLI has now converged at iteration 20 under an 8.06 GB peak RSS. The restartable six-stage workflow has not yet been rerun with the converged setting. See the [source-atlas report](docs/reports/toads-source-vandy-atlas-comparison-2026-10-03.md) and [brain-resolution convergence report](docs/reports/toads-brain-resolution-converged-2026-10-04.md).

In that bounded base run, both TOADS passes reached their two-iteration caps (terminal ratios 1.10402799 and 0.571384311); that six-stage result remains nonconverged. The source-multigrid port now matches the selected Java implementation on the full `(163,207,156)` crop: maximum component error is `6.77e-6`, below the `5e-5` gate. The base workflow recomputed core and refinement from this solver through its resume path; all stage checks passed. The compact Java/Python TOADS comparison matches startup and the first cycle closely; after 12 iterations, final label Dice is at least 0.9827 on all 11 labels, although energy ratios diverge after iteration 2. These selected synthetic comparisons do not establish whole-workflow numerical acceptance.

The source-order loader change passed its focused regression command (**14 passed in 9.64s**, code commit `4b57532`). After the multigrid port, focused GVF/plus-GVF/CLI/core/restart tests passed (**19 passed in 77.68s**), then the full suite passed on `ea31121` (**1,141 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 801.80s**). Peak sampled RSS was **15.760 GB** under the 16 GB watchdog cap, so later cycles should use targeted suites where appropriate. With source-matched padding 10, compact source-atlas TOADS converged at iteration 30 under 1.23 GB RSS, but Java/Python traces still diverge after iteration 2; see the [compact convergence report](docs/reports/toads-50-iteration-compact-2026-10-04.md) and [source comparison](docs/reports/toads-source-vandy-atlas-comparison-2026-10-03.md). The production-CLI brain-resolution TOADS run then converged at iteration 20 under 8.06 GB peak RSS; output changed in 3,536 voxels between the iteration-12 and iteration-20 Python runs, with Dice above 0.982 on every label. This is iteration-stability evidence, not Java equivalence. See the [brain-resolution convergence report](docs/reports/toads-brain-resolution-converged-2026-10-04.md). The next milestone is the restartable six-stage base workflow with both TOADS passes allowed to converge. No real-MRI validation is recorded. See the [full-suite report](docs/reports/full-pytest-2026-10-04-gvf-multigrid.md).

**Next milestone:** rerun the restartable six-stage base workflow with converged brain-resolution TOADS settings, then compare complete synthetic outputs with the Java workflow. Keep Python iteration stability, source equivalence, whole-workflow numerical acceptance, and real-data validation distinct. See [plan.md](plan.md) for row-level progress.

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
