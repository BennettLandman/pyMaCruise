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

**TOADS convergence blocker:** the registered loop now defaults to the selected wrapper's atlas range `0.2` (`atlasScale`/`priorScale`) and retains its `doGrowing=true` terminal-growth behavior. A corrected Java source harness now uses the same masked synthetic T1 and 35×33×31 working grid as Python. Its initial atlas-registration membership array matches Python exactly (SHA-256 `03fafc…`), but later membership states diverge. Python's 50-cycle diagnostic took 428.9s and stopped at the cap with energy ratio `0.375456274`; it did not meet `1e-4`. The corrected compact-atlas Java harness also fails to converge and reaches NaNs by iteration 48, so it is a state-comparison harness, not valid evidence that the source workflow converges. Earlier Java completion/convergence figures were caused by a malformed descriptor and are superseded. See the [source-runtime probe](docs/reports/toads-vandycruise-synthetic-runtime-2026-10-02.md).

**Latest full suite:** `./.venv/bin/pytest -q` → **1,130 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,213.79s** on tested commit `f3b874b`; see [the full-suite report](docs/reports/full-pytest-2026-10-02-toads-kernel-optimization.md). It includes the restartable synthetic base workflow, exact scalar-oracle parity for relation factors, Java-referenced tissue-energy outputs, and the Python-only state trace. It is not whole-workflow Java equivalence or real-data validation.

**Next integration milestone:** isolate the first divergence after the matching registration memberships, with emphasis on initial membership/statistics and the reduced synthetic atlas profile. Then make the Python synthetic base run converge within a defensible cap before advancing plus. The base workflow currently completes and resumes with manifest checks, but its two-iteration TOADS stages are deliberately marked capped and do not establish base acceptance. Synthetic integration, source-runtime behavior, source equivalence, and real-data validation remain separate evidence categories. Latest full suite: `1,130 passed, 0 failed, 0 skipped, 11 expected RAS warnings` on `f3b874b`; the subsequent atlas-range and trace updates have focused results only.

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
