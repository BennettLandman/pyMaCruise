# Registered TOADS next-loop exact update (P18cybi)

## Ported behavior

Continue a registered membership cycle that did not stop through the next
loop's carried-degree decision and exact statistics update. The adapter
consumes the degree from the preceding update, uses the current branch's
float32 energy ratio and `sqrt(maximum_distance)`, and routes through the
existing multi-channel exact field/statistics/membership implementation. It
preserves the unity or fitted field as attached state, then updates exact
centroids, deviations, image variance, centroid similarity, and memberships
on the next thinning result.

The selected adapter handles correction-enabled scheduling and the
correction-disabled degree-zero path. It refuses to continue after a selected
stopping branch. This is loop-state composition; the numerical methods remain
the already-port-tested implementations.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_first_iteration_cycle.py
```

Result: **1 passed, 0 failed, 0 skipped in 27.01 seconds**.

The committed compact RAS fixture exercises the source path through the first
exact update, rigid registration, growth/relation refresh, iteration-two
thinning and memberships, a second registration/control branch, and the
following exact update. It verifies the degree transition against the selected
gate and compares centroids, deviations, image variance/standard deviation,
similarity, memberships, denominator, and objective exactly with direct
selected-kernel calls. Positive-degree field math is covered by the separate
Java-formula tests in P18cybb/P18cyba; this test records the integration's
actual source-scheduled degree rather than assuming one.

This is synthetic composition, not a combined Java wrapper run or full TOADS
equivalence. Full-suite evidence is pending.
