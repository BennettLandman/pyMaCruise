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
  --registration-max-iterations 1 --toads-iterations 2 \
  --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1
```

Ten and twelve TOADS iterations on the 36×45×34 synthetic atlas case reduce the objective from 100.35 to about 58.3–58.6, but do not reach the `1e-4` convergence threshold. The brain-resolution integration uses two TOADS iterations per pass to validate the full route within bounded runtime; its terminal ratios are 1.10341454 and 0.570566177, with both passes reporting the iteration cap and no convergence. See the [memory investigation](docs/reports/toads-memory-bounding-2026-10-03.md), [folded-atlas integration report](docs/reports/toads-folded-vandy-atlas-base-integration-2026-10-03.md), and [earlier compact source-runtime comparison](docs/reports/toads-vandycruise-synthetic-runtime-2026-10-02.md).

**Focused validation:** on `139e238`, the 13 selected membership/registered-loop tests passed on the Python fallback. After installing the declared Numba dependency, the same Java-reference cases passed with the dispatcher exercised; ordered centroid/statistics cases passed **26 tests in 6.71s**, and multichannel registration handoff passed **20 tests in 17.86s**. After progress telemetry and queue preflight changes, registered-cycle/base restart cases passed **11 tests in 48.63s** and growth queue Java-reference cases passed **3 tests in 0.97s**. The 16 GB watchdog and its over-limit/cleanup tests are covered in the [runtime report](docs/reports/toads-memory-bounding-2026-10-03.md). The last full-suite run remains **1,130 passed, 0 failed, 0 skipped, 11 expected RAS warnings** on `f3b874b`; it predates current TOADS changes. No real-MRI validation has been recorded.

**Next integration milestone:** align TOADS startup/preprocessing grids with the Java reference, then advance the brain-resolution base workflow beyond the two-iteration smoke cap and measure convergence and stage-level numerical agreement. The atlas-support scan now uses a source-ordered Numba kernel; a full brain-resolution two-iteration TOADS probe on the updated code completed under the watchdog and was bitwise identical to the previous classification output. Its TOADS stage took 633.0 seconds and peaked at 5.941 GB sampled RSS. After TOADS behavior is aligned, proceed to plus-workflow integration. Keep synthetic integration, component source matches, whole-workflow source evidence, and real-data validation distinct.

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
