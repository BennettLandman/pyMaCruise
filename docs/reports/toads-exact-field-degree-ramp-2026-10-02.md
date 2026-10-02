# TOADS exact-field degree ramp (P18cyaz)

## Ported behavior

`toads_post_thinning_exact_field_update` now models the selected exact-field
continuation after a non-stopping thinning state. It receives the previous
degree as explicit algorithm state, compares the prior energy ratio with the
strict float32 threshold `sqrt(maxDistance)`, advances at most one degree,
fits the newly selected degree 1–3 exact field, and preserves the attached
field unchanged when the selected degree remains zero. The existing degree-one
adapter continues to reject branches that do not schedule degree one.

This follows the selected short-constructor path. It does not add degree 4,
spline correction, or the rest of the centroid/membership/registration loop.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_post_thinning_exact_field_ramp.py tests/test_toads_post_thinning_exact_field.py tests/test_toads_exact_polynomial_field_degrees23.py tests/test_toads_post_thinning_exact_iteration.py
```

Result: **15 passed, 0 failed, 0 skipped in 43.29 seconds** on implementation
and test commit `11e326613f2c777a82033d3df17f79118ec06ddf`.

At a below-threshold ratio, current degrees 0, 1, and 2 schedule 1, 2, and 3.
Each completed field is compared with the selected Java formula harness on
the committed 21³ folded-cortex RAS phantom; degree-2/3 maximum field error is
bounded by `2e-6`. The degree-zero no-op leaves its input field unchanged.
Existing degree-one post-thinning and carried-iteration tests remain green.

The full repository suite has not yet been rerun on this commit, so P18cyaz
remains Started.
