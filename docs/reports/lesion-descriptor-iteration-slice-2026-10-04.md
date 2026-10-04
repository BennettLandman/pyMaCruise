# LesionTOADS descriptor data-to-membership slice

This bounded composition connects the selected active 13-class raw energy
dispatch with the descriptor-aware, source-ordered membership updater. It
consumes already estimated centroids and channel weights, fixed shape and
similarity maps, and pre-resolved direct/prior relation choices. It does not
estimate statistics, rebuild relation maps, apply lesion-distance modulation,
refit correction fields, or control convergence.

The pinned Java harness uses the bundled 13-class order, two-channel
synthetic intensities, and the active source switch cases. With unit shapes,
direct relations, zero smoothing, and zero prior coefficient, raw class
energy, all 13 normalized class memberships, shared lesion energy and
membership, and the denominator match their Java float32 bit patterns. This is
a fixed-state single-voxel source-formula match, not a connected iteration or
workflow-equivalence result.

Focused command:

```text
./.venv/bin/pytest -q tests/test_lesion_descriptor_iteration_slice.py tests/test_lesion_active_class_energy.py tests/test_lesion_descriptor_membership.py tests/test_lesion_objective_membership_pass.py
```

Result: **9 passed, 0 failed, 0 skipped in 0.62 seconds**.

Full suite on tested commit `3aadd8b18e3d80fdb90a5f5fe2ed9a9faa0bf78d`:
**1,150 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 846.38s**,
with peak sampled RSS **12.444 GB** under the 24 GB watchdog. See the
[full-suite report](full-pytest-2026-10-04-descriptor-iteration-slice.md).

Statistics refresh, iterative stopping, complete spatial/atlas integration,
full workflow comparison, and real-MRI validation remain open.
