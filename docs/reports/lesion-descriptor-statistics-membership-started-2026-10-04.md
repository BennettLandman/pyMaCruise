# Descriptor statistics-to-membership transition (started)

This adapter now runs the selected post-thinning order for the bundled
13-class descriptor: class/lesion/black-hole centroid refresh, variances,
channel weights, centroid similarity, active class raw terms, membership
updates, then the energy ratio. It maps descriptor `obj` to the statistics
module's `object` token, while retaining original descriptor topology for the
active energy/membership dispatch.

All model state is still synthetic or caller-supplied: hard segmentation,
initial memberships, prior centroids, and shape/relation maps. This does not
run thinning, rebuild relations, update correction fields, apply lesion
distance modulation, or continue the outer iteration. `use_lesion_weight=True`
is explicitly rejected because the descriptor raw-energy adapter has not yet
implemented the selected FLAIR lesion PowerTable branch.

Focused command:

```text
./.venv/bin/pytest -q tests/test_lesion_descriptor_statistics_membership.py tests/test_lesion_statistics_membership_bridge.py tests/test_lesion_statistics_refresh.py tests/test_lesion_descriptor_iteration_slice.py
```

Result: **6 passed, 0 failed, 0 skipped in 0.68 seconds**. The new 5×5×5
synthetic test compares this adapter against the same already-validated
statistics, similarity, and descriptor transition components composed
independently, and confirms inputs are not mutated. It is integration
plumbing evidence only, not a new Java harness match. Next, compare refreshed
intermediates to Java on a 13-class fixture and implement the active FLAIR
lesion-weight option before expanding the supported wrapper profile.
