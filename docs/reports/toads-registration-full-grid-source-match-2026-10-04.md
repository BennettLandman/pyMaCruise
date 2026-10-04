# TOADS full-grid registration source match — 2026-10-04

## Scope

This report corrects the earlier diagnosis in
[the two-iteration boundary report](toads-vandy-iteration-boundary-2026-10-04.md).
The Python post-refresh adapter had limited exact registration coefficients
to per-class support bounds. The selected Java
`MultipleToadDeformableAtlas.computeExactRegistrationCoefficients` loops across
the complete image grid for every registered shape; it does not apply those
bounds. The adapter now supplies full-grid bounds for coefficient calculation
and still computes shape bounds after registration for transformed-output
consumers.

## Java source evidence

`tests/java_reference/ToadsVandyPostRefreshRegistrationHarness.java` configures
the selected Java atlas class and captures the registration input, normalized
gradient/Hessian, baseline energy, and final transform from the committed
synthetic Vandy arrays. Regenerate its reference with:

```bash
PYTHONPATH=src ./.venv/bin/python tests/java_reference/generate_toads_vandy_post_refresh_reference.py
```

At the post-refresh starting transform and Java iteration-1 memberships:

- The Python normalized gradient and Hessian each match Java bit-for-bit.
- The Python baseline energy differs by one float32 ULP.
- Running the Python registration update on those same Java inputs returns
  the exact Java transform bits.
- The integrated Python post-refresh transform is within one ULP of that
  source output; the first registered-loop transform then matches Java
  bit-for-bit.

The prior claim that a `1.93e-8` post-refresh transform difference caused the
first loop split is superseded. That difference came from the support-bounded
coefficient scan and is corrected in implementation commit `917c15d`.

## Remaining synthetic differences

The transform correction does not eliminate all iteration drift:

| State | Python vs selected Java |
|---|---:|
| First updated memberships | mean absolute error `1.2715436e-10`; maximum `4.529953e-6`; 16 values exceed `1e-6` |
| Registered-loop transform | Bit-identical |
| Post-growth labels | 55 of 162,180 voxels differ |
| Closest-label state | 843 of 648,720 entries differ |
| Next thinning | 42 of 162,180 voxels differ; minimum class Dice `0.99558791` |
| Iteration-2 centroids | MAE `7.179651e-8`; maximum `7.152557e-7` |
| Iteration-2 standard deviations | MAE `1.8164747e-5`; maximum `0.000109164044`; five classes exceed `1e-6` |
| Iteration-2 memberships | MAE `5.3438457e-6`; maximum `0.39456701`; 2,006 values exceed `1e-5` |
| Iteration-2 objective | Python `43.74055862`; Java `43.742954`; absolute difference `0.00239538` |

The first remaining discrete difference appears after the now-matching
registration handoff, in post-growth labels and closest-label state. Iteration
2 remains visibly divergent locally; no membership or whole-TOADS equivalence
is claimed. No clinical image was used.

## Validation

Focused registration group on implementation commit `917c15d`:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q tests/test_toads_registration_sampling.py tests/test_toads_vandy_iteration_boundary.py tests/test_toads_post_registration_template_refresh.py tests/test_toads_post_registration_branch.py tests/test_toads_alignment_startup.py
```

Result: **32 passed, 0 failed, 0 skipped in 12.66s**.

Full suite on the same implementation commit, under the 24 GB RSS watchdog:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,159 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
871.73s**. Peak sampled RSS was **14.655 GB**; the watchdog did not terminate
the run. This is synthetic selected-source slice evidence, not whole-workflow
equivalence or real-data validation.
