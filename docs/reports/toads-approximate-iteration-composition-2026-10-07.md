# Alternate TOADS approximate iteration composition

## Scope

Added `run_default_toads_approximate_iteration` to connect the alternate
default-atlas `approx=4` path after a supplied thinning state. It calls the
source-ordered approximate centroid update, computes assigned-class deviations
using the refreshed centroids, then calls the four-slot approximate membership
pass with those refreshed centroids and deviations. The adapter returns each
stage's accumulators and output state and does not mutate its input arrays.

This is not the short `MedicAlgorithmToads08` production path: that constructor
uses exact memberships (`useApprox=false`). It is also not a Java-composite
iteration comparison or a complete approximate TOADS loop. Centroid,
deviation, and membership methods have their own selected-Java reference
fixtures; this new composition test compares the adapter's outputs
bit-for-bit with direct calls to those already-tested methods.

## Validation

Command:

```sh
./.venv/bin/pytest -q --tb=short \
  tests/test_toads_default_approximate_iteration.py \
  tests/test_toads_approximate_centroids.py \
  tests/test_toads_approximate_variances.py \
  tests/test_toads_default_approximate_membership_pass.py
```

Result: **7 passed, 0 failed, 0 skipped in 19.27s** on implementation commit
`fce1e9021983e4959aff18b7a5affc32e6e441f5`. This was a focused component and
composition run; the full suite was not rerun for this adapter-only addition.
No selected Java source lines were newly represented in this cycle.

## Next integration step

Carry the composed approximate membership state into the source-ordered
registration/control and next thinning/growing transition, then compare a
complete alternate-constructor iteration boundary to an instrumented Java
reference. Keep this evidence separate from the selected production wrapper
and the brain-resolution base workflow comparisons.
