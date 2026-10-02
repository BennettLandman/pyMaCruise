# Registered TOADS iteration-two exact update (P18cybg)

## Ported behavior

Continue the P18cybf registration/growth/relation/next-thinning result through
the next source-ordered field decision and exact-statistics update. The
wrapper's energy ratio from iteration one is compared with
`sqrt(maximum_distance)` to advance the carried degree when inhomogeneity
correction is enabled. Degree zero leaves the attached unity field unchanged;
the code then recomputes exact centroids, class deviations, image variance,
centroid similarity, and 11-class memberships on the iteration-two thinned
state. The correction-disabled branch keeps degree zero and skips the field
stage, while running the same statistics and membership update.

The implementation routes through the existing multi-channel carried-degree
adapter. It adds no new numerical kernel. The wrapper settings (maximum
distance, exponents, modality and variance weights, shape and smoothing
priors) remain explicit inputs.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_first_iteration_cycle.py
```

Result: **1 passed, 0 failed, 0 skipped in 16.88 seconds**.

The committed compact RAS phantom traverses actual initial registration,
iteration-one exact memberships and rigid registration, continuing growth,
relation refresh, shape transformation, iteration-two thinning, and the
iteration-two exact update. Direct calls to the existing centroid, class
variance, image variance, centroid similarity, and membership kernels match
all resulting arrays, denominators, and objective exactly. The test also
checks the correction-disabled branch and unchanged unity field.

With `maximum_distance=0.04`, the measured iteration-one energy ratio does
not fall below the degree threshold `sqrt(0.04)=0.2`; the correct result here
is degree zero with no refit. Positive-degree fitting remains validated in
the carried-degree Java formula comparisons (P18cybb/P18cyba), but this
integration slice does not claim an end-to-end positive-degree refit.

This is synthetic composition, not a combined Java wrapper execution or full
TOADS equivalence. Full-suite evidence is pending.
