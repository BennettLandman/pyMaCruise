# LesionTOADS descriptor data-to-membership slice (started)

This slice composes P39zzzzzzbd's active raw data-term dispatch with
P39zzzzzzbe's descriptor-aware source-ordered membership updater. It consumes
already estimated centroids and channel weights, fixed shape/similarity maps,
and direct/prior relation decisions. It does not estimate statistics, rebuild
relations, apply lesion distance modulation, refit correction fields, or
control convergence.

The pinned Java harness uses the bundled 13-class order, two-channel synthetic
intensities, fixed centroids, and the selected active switch cases. With unit
shapes, direct relations, zero smoothing, and zero prior coefficient, raw
class energy, all 13 normalized class memberships, shared lesion energy and
membership, and the membership denominator match their float32 Java bit
patterns.

Focused command:

```text
./.venv/bin/pytest -q tests/test_lesion_descriptor_iteration_slice.py tests/test_lesion_active_class_energy.py tests/test_lesion_descriptor_membership.py tests/test_lesion_objective_membership_pass.py
```

Result: **9 passed, 0 failed, 0 skipped in 0.62 seconds**. This is a bounded
single-voxel source-formula match with fixed statistical state, not an
iterative LesionTOADS or end-to-end workflow result. The full suite remains
pending.
