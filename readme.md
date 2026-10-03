# PyMaCRUISE

Python-native implementation of the MaCRUISE cortical segmentation workflow. The project reuses established imaging libraries and ports the source-specific algorithms and orchestration that are not available as suitable Python components.

[plan.md](plan.md) is the detailed, step-by-step port and validation plan. This README keeps only the current project status and stable project guidance; historical validation reports remain under [docs/reports](docs/reports/).

## Current status — 2026-10-02

**The restartable synthetic base run completes every stage, but TOADS convergence remains open.** `pymacruise base-run` composes ANTsPy N4 → registered TOADS → Dura → second TOADS → MaSeg cleanup → ACE/GVF/TGDM → MaSeg refinement. The command below wrote `inner.nii.gz`, `central.nii.gz`, `outer.nii.gz`, `SEG_refine.nii.gz`, the rule summary, stage logs, and `run.json`; all six stages passed their checks. A `--resume` rerun verified all six and skipped recomputation. The manifest keeps both two-iteration TOADS stages marked capped and `converged=false`.

**Runnable synthetic command** (from repository root; expected inputs are committed fixtures):

```bash
./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/toads_registered_first_iteration_cycle/t1.nii.gz \
  --brain-mask tests/data/synthetic/toads_registered_first_iteration_cycle/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/toads_registered_field_refit/atlas-11.txt \
  --output-dir build/synthetic-base-run \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --toads-iterations 2 \
  --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1
```

The runner's output is a synthetic integration artifact, not an accepted segmentation or whole-workflow equivalence claim. See the [restartable base-run report](docs/reports/restartable-base-workflow-synthetic-2026-10-02.md).

**Newly source-compared TOADS behavior:** membership uses the fixed atlas Intensity Variance Factor while refreshed class standard deviations feed centroid similarity. A two-update synthetic Java harness matches both objective float32 bit patterns and the source energy-ratio bits. This is selected update-sequence equivalence only. Focused results: 21 TOADS integration tests passed in 74.52s on `224572f`; two-update parity test passed in 13.92s on `029255a`; base restart test passed in 80.02s on `224572f`.

**Latest full suite:** `./.venv/bin/pytest -q` → **1,125 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,265.88s** on tested commit `029255a`; see [the full-suite report](docs/reports/full-pytest-2026-10-02-toads-variance-factor-correction.md). These are synthetic regression and bounded source-comparison results, not whole-workflow Java equivalence or real-data validation.

**Remaining blocker and next milestone:** with the factor correction, the 13³ loop still alternates and misses `1e-4` after 40 iterations (56.712s; terminal ratio 0.0397632). Shape-prior and intensity-scale probes also failed to converge within 10 iterations. The next TOADS milestone is a source comparison of the complete outer-loop state handoff through thinning, growing, registration, and the next membership update; address the alternating behavior and accelerate the ordered membership sweep before raising workflow caps. Then advance plus with synthetic ANTsPy rigid FLAIR-to-T1 registration and the LesionTOADS handoff. The plus downstream test currently starts from crop-ready TOADS products; full LesionTOADS and real-MRI validation remain open. See [the TOADS variance-factor and convergence report](docs/reports/toads-variance-factor-and-convergence-2026-10-02.md) and [the plus integration report](docs/reports/plus-downstream-synthetic-2026-10-02.md).

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
