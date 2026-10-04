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

**Validation:** the latest focused run, `.venv/bin/pytest -q tests/test_toads_loop_cli.py tests/test_base_workflow_restartable.py`, passed **7 tests in 68.72s**. The last full suite passed **1,138 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 793.64s** on `ad34321`, with peak sampled RSS **13.960 GB** under the 16 GB watchdog; it has not been rerun since. See the [full-suite report](docs/reports/full-pytest-2026-10-03-brain-resolution.md) and [runtime and Java-comparison report](docs/reports/toads-memory-bounding-2026-10-03.md). No real-MRI validation has been recorded.

**TOADS convergence and source comparison:** on the compact `(36,45,34)` atlas-derived synthetic grid, Python reached its `1e-4` energy-ratio threshold at iteration 12 (`7.53300919e-5`) in 38.8 seconds, with 0.593 GB sampled peak RSS. The source harness uses one post-refresh rigid update (`alignIter=1`); Python had been compared at its five-update default, so the earlier transform mismatch was an invalid comparison. With the source setting, transforms agree within about `5e-7` and transformed shape priors within `7.4e-6`. The first captured source/Python divergence is after the first loop's thinning/statistics/membership update: 23,240 voxels differ in at least one membership by more than `1e-6`; energy ratios are `0.60919` vs `0.65671`. Final ordinal-label voxel agreement is 98.06%, but only 2 of 11 per-label Dice scores meet 0.97. This remains a source mismatch, not equivalence. At brain resolution `(172,216,164)`, a 15-iteration continuation remained nonconverged at `0.0139147937`; it took 2,103 seconds and peaked at 7.995 GB sampled RSS under the 16 GB watchdog. Ratios were oscillatory (0.02161 at iteration 10, 0.23222 at 11, 0.01391 at 15), so compact-grid convergence does not transfer to the larger grid. The final brain-resolution uint8 classification is RAS with the expected 11 labels and SHA-256 `8673a38edec9ee25b3b5ec75d91333330e48f6084a6ad893043422173f36ecee`. No real-MRI validation is claimed. See the [runtime and Java-comparison report](docs/reports/toads-memory-bounding-2026-10-03.md).

A newer source-atlas diagnostic compared the original VandyCRUISE raw atlas directly with its bundled NIfTI conversion; 35 tested atlas payloads were voxel-exact after source-order or source-header restoration. Java and Python still crop to different x extents by one voxel. On their crop-derived shared region, transformed shape priors differ by mean absolute `0.00986` (maximum `1.0`), startup memberships by `0.00737`, and transformed topology agrees at 92.04% of voxels. This is synthetic source-comparison evidence of an unresolved mismatch, not equivalence. Details are in the [source-atlas diagnostic](docs/reports/toads-source-vandy-atlas-comparison-2026-10-03.md) and the [earlier provisional-atlas diagnostic](docs/reports/toads-first-loop-source-diagnostic-2026-10-03.md).

The atlas conversion reader now restores source voxel order when RAS normalization permutes axes, as well as when it only flips axes. Its focused atlas suite passed **7 tests in 0.84s**. This covers the general conversion API; it does not certify the temporary diagnostic atlas conversion.

**Next integration milestone:** align the Java/Python startup crop and compare constructor labels, centroids, variances, and initial memberships on the same voxel grid. Then continue first-loop comparisons before longer convergence cases. Keep synthetic integration, component source matches, whole-workflow source evidence, and real-data validation distinct.

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
