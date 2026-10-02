# First TOADS exact loop iteration from unity-field startup (P18cybe)

## Ported behavior

Continue the source-faithful registered startup after its initial membership
pass and first thinning. The wrapper's initial `Eratio=1.0` and degree zero
leave the field unchanged during the first loop's correction stage. The code
then updates exact centroids from the thinned memberships/labels, refreshes
class deviations and image-wide variance, recomputes centroid similarity, and
runs the next exact membership pass using the registered shape and relation
maps.

The field remains attached by identity and is still all ones after this
degree-zero stage. This adapter stops before atlas registration, growth, and
outer-loop control for the new membership output.

## Validation

Focused command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_field_refit.py::test_registered_unity_startup_continues_through_first_degree_zero_exact_loop
```

Result: **1 passed, 0 failed, 0 skipped in 25.32 seconds**.

The committed RAS synthetic atlas/T1/mask case runs through registration,
initial unity-field membership, first thinning, and the degree-zero exact
loop update. The test verifies the carried field is unchanged, degree stays
zero, stage order matches the wrapper, and the new membership maps have the
same default 11-class grid. It also recomputes every stage directly with the
existing selected kernels and compares centroids, class deviations, image
variance and standard deviation, similarity, membership maps, denominators,
and objective exactly. Those kernels also have separate selected Java-reference
tests. This is synthetic composition, not a combined Java wrapper run or full
TOADS equivalence.
