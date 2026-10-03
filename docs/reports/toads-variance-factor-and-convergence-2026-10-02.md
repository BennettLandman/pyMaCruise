# TOADS atlas variance-factor correction and convergence diagnostic

## Source behavior corrected

VandyCRUISE `MultipleToadSegmentation` keeps atlas `varweight` separate from
class `stddev`. The selected `AlgorithmMultipleToads` loop refreshes class
standard deviations for centroid similarity, but the selected membership
calculation uses the atlas Intensity Variance Factor. The call that would
re-estimate variance weights is commented out in the selected outer loop.
Python had passed measured deviations into membership at several composed
entry points. The descriptor reader and those integrations now preserve the
fixed atlas factors through each exact membership update.

## Validation

- `./.venv/bin/pytest -q tests/test_toads_intensity_variance_factors.py tests/test_toads_initial_statistics_membership.py tests/test_toads_initial_statistics_thinning.py tests/test_toads_initial_thinning_statistics_membership.py tests/test_toads_registered_initial_membership.py tests/test_toads_registered_atlas_membership_integration.py tests/test_toads_registered_first_iteration_cycle.py tests/test_toads_registered_field_refit.py tests/test_toads_initial_registration_control.py` — **21 passed** in 74.52s on implementation commit `224572f`.
- `./.venv/bin/pytest -q tests/test_toads_two_exact_updates_source_order.py` — **1 passed** in 13.92s on commit `029255a`. The Java harness and Python composition match both exact objective float32 bit patterns and the energy-ratio bits; it separately checks fixed factor bits and refreshed deviation bits.
- `./.venv/bin/pytest -q tests/test_base_workflow_restartable.py` — **1 passed** in 80.02s on commit `224572f`.
- Full repository suite — **1,125 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,265.88s** on commit `029255a`; see the [full-suite report](full-pytest-2026-10-02-toads-variance-factor-correction.md).

A manual `pymacruise base-run` invocation on `224572f` completed all six stages. It wrote the final native-grid RAS `.nii.gz` products (`inner`, `central`, `outer`, and `SEG_refine`), the refinement rule summary, stage intermediates, logs, and `run.json`. The manifest recorded all six stages complete and both two-iteration TOADS calls as iteration-capped with `converged=false`. A rerun with `--resume` verified each stage and skipped all six.

## Remaining convergence and runtime blocker

The 13³ synthetic registered-loop trace was rerun with the corrected descriptor
factors:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --output build/variance-corrected-toads-loop-40.nii.gz \
  --maximum-iterations 40 --convergence-threshold 1e-4
```

It took 56.712s and stopped at the 40-iteration cap with ratio 0.0397632,
not at the 1e-4 threshold. The exact objective alternates between two bands;
this is not explained by intensity scale alone. Ten-iteration probes with
`descriptor`, `uniform-one`, and `topology-one-hot` shapes, plus `original-scaled`
intensities, all missed 1e-4. Descriptor and scaled-intensity controls emitted
the same terminal label hash. These changes were diagnostic perturbations,
not source-faithful alternatives or claims about Java convergence.

The Java/Python two-update comparison validates the selected membership and
energy-ratio seam, not the full outer loop with thinning, growing, registration,
and per-iteration state transfer. That whole-loop state transition is the
next TOADS integration target. The source-ordered membership sweep remains a
runtime hotspot on larger volumes; no JIT dependency has been added. Real-MRI
validation remains deferred.
