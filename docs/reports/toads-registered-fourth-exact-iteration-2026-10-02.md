# TOADS iteration-four carried-degree exact update (P18cybk)

## Ported behavior

When the P18cybj registered membership branch continues, pass its next-thinning
state through the source energy-ratio degree decision and the following exact
statistics/membership update. Preserve attached correction fields and compare
centroids, class deviations, image variance, centroid similarity, exact
memberships, denominator, and objective against direct selected-kernel calls.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_first_iteration_cycle.py
```

Result: **1 passed, 0 failed, 0 skipped in 31.16 seconds**.

The test follows the continuing branch conditionally. It checks the scheduled
degree and field-refit/no-op behavior, confirms the attached fields, then
recomputes every update using the direct selected kernels and compares outputs
exactly. This extends the same compact synthetic RAS workflow; it is not a
combined Java wrapper execution or full TOADS equivalence. Full-suite evidence
is pending, so P18cybk remains Started.
