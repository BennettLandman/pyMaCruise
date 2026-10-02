# TOADS initial class variances (P18cyam)

**Implementation commit:** `1aaa5419dc4a97086f6f80cd3e1ceb18603b58a2`  
**Validation date:** 2026-10-02

The atlas-backed startup composition now computes the selected exact
assigned-class deviations after initialization of memberships and relation
state. The selected MASK and optimized CSF/GM/WM branches share the source
formula: scan the full image grid in x/y/z order, include only voxels assigned
to that class, accumulate `(membership² + 1e-30f) * residual²`, and return
`sqrt(numerator / denominator)` or zero when no voxels are assigned. The
initial correction fields contain only ones, so omitting an explicit unity
field array preserves the operation without allocating another full volume.

Focused command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_initial_relation_state.py tests/test_toads_relation_list.py tests/test_toads_boundary_distances.py tests/test_toads_closest_labels.py tests/test_toads_relation_state_refresh.py tests/test_toads_optimized_class_variance.py
```

Result: **38 passed, 0 failed, 0 skipped in 7.78 seconds** on the implementation
commit above. The synthetic atlas composition verifies the variance stage and
output shape after startup; the underlying formula is checked against the
selected Java reference, including nonunit field and empty-class cases. This
does not yet compute image variance, centroid similarity, or the initial
membership/energy sweep.
