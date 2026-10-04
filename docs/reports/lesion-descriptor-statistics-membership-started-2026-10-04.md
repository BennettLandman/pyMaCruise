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
distance modulation, or continue the outer iteration. The optional
`use_lesion_weight=True` branch routes the descriptor's FLAIR channel and
per-channel lesion PowerTable exponent through statistics refresh and raw
LESION energy calculation.

Focused command:

```text
./.venv/bin/pytest -q tests/test_lesion_active_class_energy.py tests/test_lesion_descriptor_statistics_membership.py tests/test_lesion_data_energy_maps_special.py
```

Result: **8 passed, 0 failed, 0 skipped in 0.69 seconds**. The new 5×5×5
synthetic test compares this adapter against the same already-validated
statistics, similarity, and descriptor transition components composed
independently, checks weighted and unweighted FLAIR lesion-energy modes, and
confirms inputs are not mutated. The raw lesion term reuses the previously
source-matched scalar PowerTable helper; this orchestration test itself is
integration plumbing evidence, not a new Java harness match. Missing exponents
and out-of-range FLAIR channels fail explicitly. Next, compare refreshed
intermediates to Java on a 13-class fixture and validate the complete wrapper
profile. The row remains Started pending a full-suite evaluation.
