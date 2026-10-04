# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-03

**The synthetic base MaCRUISE route now completes and resumes at brain resolution** on a deterministic `(172, 216, 164)` atlas-derived phantom. All six stage checks pass; native-grid final level sets and refined labels are RAS; the manifest records all stage hashes and explicitly records non-convergence. A second invocation verifies/skips all six stages. The run peaked at 6.615 GB sampled RSS under the 16 GB watchdog; resume verification peaked at 0.123 GB. This is synthetic integration evidence, not a converged or source-equivalent segmentation and not real-data validation. Full run evidence and the reproducible command are in the [brain-resolution runtime report](docs/reports/toads-memory-bounding-2026-10-03.md).

**Runnable synthetic command** (from repository root):

```bash
./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/toads_downsampled_vandy_atlas/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_downsampled_vandy_atlas/brain-mask.nii.gz \
  --maseg tests/data/synthetic/toads_downsampled_vandy_atlas/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_downsampled_vandy_atlas/atlas.txt \
  --output-dir build/downsampled-vandy-base-run \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 \
  --registration-max-iterations 1 --toads-iterations 2 \
  --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1
```

Ten and twelve TOADS iterations on the 36×45×34 synthetic atlas case reduce the objective from 100.35 to about 58.3–58.6, but do not reach the `1e-4` convergence threshold. The brain-resolution integration uses two TOADS iterations per pass to validate the full route within bounded runtime; its terminal ratios are 1.10341454 and 0.570566177, with both passes reporting the iteration cap and no convergence. See the [memory investigation](docs/reports/toads-memory-bounding-2026-10-03.md), [folded-atlas integration report](docs/reports/toads-folded-vandy-atlas-base-integration-2026-10-03.md), and [earlier compact source-runtime comparison](docs/reports/toads-vandycruise-synthetic-runtime-2026-10-02.md).

**Validation:** after the source-order atlas loader change, the focused atlas/TOADS command passed **14 tests in 9.64s** on `4b57532`. The last full suite passed **1,138 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 793.64s** on `ad34321`, with peak sampled RSS **13.960 GB** under the 16 GB watchdog; it has not been rerun since. See the [full-suite report](docs/reports/full-pytest-2026-10-03-brain-resolution.md) and [runtime and Java-comparison report](docs/reports/toads-memory-bounding-2026-10-03.md). No real-MRI validation has been recorded.

**Current compact source comparison:** a 12-iteration Java/Python run on the `(36,45,34)` Vandy-atlas synthetic case gives 99.82% final voxel agreement and Dice ≥0.9827 for all 11 labels. Neither loop converged: Python's terminal energy ratio was `0.0003041`, Java's was `0.02334`; their ratios begin to diverge after iteration 2. The low-memory Python run completed in 37.7 seconds with history and per-event capture disabled. Retaining those diagnostic events caused an earlier 12-iteration trace to exit with code 137, so the longer run used the low-memory setting. Brain-resolution convergence has not been rerun after the source-order atlas fix; the previous 15-iteration `(172,216,164)` baseline ended at `0.0139148` in 2,103 seconds and 7.995 GB sampled peak RSS. No real-MRI validation is claimed. See the [runtime report](docs/reports/toads-memory-bounding-2026-10-03.md) and [source-atlas comparison](docs/reports/toads-source-vandy-atlas-comparison-2026-10-03.md).

A source-order atlas-loading fix now matches the Java/Python synthetic startup and first TOADS cycle: crop bounds and constructor topology match exactly, the first thinning labels match all 162,180 voxels, and only 67 of 1,783,980 first-update membership values differ by more than `1e-6` (maximum absolute difference `3.78e-6`). This is a selected one-cycle source-reference result; full TOADS convergence and whole-workflow numerical acceptance remain open. Details are in the [source-atlas diagnostic](docs/reports/toads-source-vandy-atlas-comparison-2026-10-03.md) and the [earlier provisional-atlas diagnostic](docs/reports/toads-first-loop-source-diagnostic-2026-10-03.md).

The atlas conversion reader now restores source voxel order when RAS normalization permutes axes, as well as when it only flips axes. Its focused atlas suite passed **7 tests in 0.84s**. This covers the general conversion API; it does not certify the temporary diagnostic atlas conversion.

**Next integration milestone:** continue source comparison after the matched first TOADS cycle, then rerun the restartable base workflow with source-ordered atlas loading. Keep synthetic integration, component source matches, whole-workflow source evidence, and real-data validation distinct.

Test totals and checklist progress are evidence bookkeeping, not measures of workflow completion. See [plan.md](plan.md) for the row-level validation status.
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
