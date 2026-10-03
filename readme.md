# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-03

**The synthetic base MaCRUISE route completes and resumes on a folded phantom derived from the bundled VandyCRUISE atlas.** All six stage checks pass, all generated NIfTI outputs are RAS, and a second invocation verifies/skips all stages. The manifest correctly says the route is not converged because both TOADS passes were capped at two iterations. This is synthetic integration evidence, not whole-workflow Java equivalence or real-data validation.

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

Ten and twelve TOADS iterations on the 36×45×34 synthetic atlas case reduce the objective from 100.35 to about 58.3–58.6, but do not reach the `1e-4` convergence threshold. The Vandy-atlas base command above now completes all six stages on that synthetic grid and verifies every stage on `--resume`; its two TOADS passes are intentionally capped at two iterations and report non-convergence. All final NIfTI outputs are RAS. The whole-brain-sized `(172, 216, 164)` synthetic workflow has been attempted three times; each reached `toads-first`, sampled RSS peaked at 4.014 GB, and no TOADS stage output was committed to the manifest before interruption. This is evidence about observed memory during the incomplete attempts, not a completed-run peak. The TOADS thinning initializer is now vectorized with source-ordered heap insertion; a fresh monitored whole-brain run remains needed to assess runtime and peak memory. See the [memory investigation](docs/reports/toads-memory-bounding-2026-10-03.md), [folded-atlas integration report](docs/reports/toads-folded-vandy-atlas-base-integration-2026-10-03.md), and [earlier compact source-runtime comparison](docs/reports/toads-vandycruise-synthetic-runtime-2026-10-02.md).

**Focused validation:** bounded loop code on `da691f5` and thinning-history follow-up on `b897a45` passed their recorded source-reference/restart checks. The current vectorized initializer and adjacent registration/base workflow tests passed **21 tests in 175.85s**, followed by **2 passed in 63.98s** for the base restart test after adding live-log assertions; a separate edge-case comparison against literal Java-order loops passed (**1 passed in 0.35s**). The complete small Vandy-atlas synthetic run produced a manifest with six completed stages; rerunning the command reported `Resume: verified` for all six. Its watchdog sampled a 0.598 GB peak. The whole-resolution fixture generator in `0708e76` produced a `(172, 216, 164)` RAS phantom. The RSS watchdog in `3c0b132`, with cleanup fix `462f032`, passed over-limit and interruption-cleanup checks. The last full-suite run remains **1,130 passed, 0 failed, 0 skipped, 11 expected RAS warnings** on `f3b874b`; it predates current TOADS changes. No complete brain-resolution Python run or real MRI validation has been recorded.

**Next integration milestone:** restart the whole-brain synthetic base run with the vectorized initializer and live TOADS phase logging, then record full-run peak RSS, stage outputs, and convergence telemetry. The watchdog samples the Python process's RSS every 0.05 seconds and terminates its process group at 16 GB; it is not a kernel-enforced cap and can overshoot between samples. Keep synthetic integration, component source matches, whole-workflow source evidence, and real-data validation distinct.

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
