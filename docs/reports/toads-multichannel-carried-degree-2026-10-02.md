# Multi-channel carried-degree TOADS iteration (P18cybb)

## Ported behavior

Added a channel-last continuation for the selected 11-class exact TOADS
iteration. It carries the previous degree, fits each channel independently in
source channel order at degrees 1–3, applies the source degree-zero field
no-op, and then refreshes centroids, class deviations, image-wide variance,
centroid similarity, and exact memberships. Positive degree updates mutate
the caller's attached field array, matching the Java wrapper's in-place field
state. Existing channel and statistics kernels are reused.

## Validation

Focused command:

```text
./.venv/bin/pytest -q tests/test_toads_multichannel_carried_degree_iteration.py
```

Result: **2 passed, 0 failed, 0 skipped in 14.48 seconds**.

The two-channel degree-three field maps agree with separate Java formula
harness runs within `2e-6`, and their returned channel order is `(0, 1)`.
Centroids, deviations, image variance, image standard deviations, similarity,
membership maps, denominators, and objective exactly match direct composition
of the existing selected Python stage kernels. A degree-zero case verifies
that both attached fields remain unchanged while statistics and memberships
continue. The folded input is committed RAS `.nii.gz` synthetic data. The
field and stage kernels have separate selected-Java reference checks; this
adapter test is not a combined Java wrapper execution.

## Scope limits

This composes a fixed-input multi-channel iteration only. Atlas transforms,
shape/relation refresh, topology decisions, growth, and outer-loop termination
remain explicit upstream/downstream stages. Full-suite evidence remains
pending.
