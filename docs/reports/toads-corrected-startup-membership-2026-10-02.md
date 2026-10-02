# Field-corrected registered TOADS startup membership (P18cyaw)

## Ported behavior

`compute_registered_toads_corrected_initial_membership` now composes the
selected startup path through the degree-three exact correction fit and the
field-aware initial statistics/exact-membership pass:

1. Fit one source-ordered correction field from each registered RAS image,
   using the initialized exact memberships, constructor labels, and centroids.
2. Recompute assigned-class deviations, image-wide variance weights, and
   centroid similarities with those fitted fields.
3. Run the selected exact-membership/data-energy update using the registered
   shape and relation state plus caller-supplied atlas/smoothing configuration.

The registered synthetic integration currently uses one T1 channel. Wrapper
priors and coefficients remain explicit inputs. This does not continue through
TOADS thinning/growing or establish a full wrapper execution.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_field_refit.py tests/test_toads_registered_atlas_membership_integration.py tests/test_toads_initial_statistics_membership.py tests/test_toads_default_atlas_dispatch.py tests/test_toads_optimized_class_variance.py tests/test_toads_optimized_variance.py tests/test_toads_object_image_variance.py tests/test_toads_optimized_similarity.py tests/test_toads_exact_polynomial_field_degrees23.py
```

Result: **35 passed, 0 failed, 0 skipped in 47.65 seconds** on implementation
and test commit `0586a7772fb3a5b56da31167f869af4c1513d0d3`.

The degree-three fit's normal matrix and right-hand side match the selected
Java formula harness exactly; its complete field is within `2e-6`, and values
outside the correction mask remain one. The integrated field-aware membership
and class-energy arrays match a direct invocation with the same run-derived
state, and memberships differ from the uncorrected control. The statistics
and membership kernels have their own selected Java-reference tests; this
adapter comparison is not an independent combined Java wrapper oracle.

The full repository suite has not yet been rerun on this commit, so P18cyaw
remains Started.
