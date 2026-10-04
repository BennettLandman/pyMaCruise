# LesionTOADS descriptor membership entry

This slice maps the bundled descriptor order to membership-update categories:
`obj` to OBJECT, `mask` to MASK, `optimized_csf` to the ordinary class branch,
and `lesion` to LESION. Optimized-CSF's raw data term is still computed by the
source-active fallback in P39zzzzzzbd. This mapping preserves separate raw
energy dispatch and membership-update dispatch.

The selected membership updater previously required one BLACKHOLE and one
OUTMASK class. Neither appears in the bundled 13-class LesionTOADS descriptor.
The kernel now permits those optional families to be absent, while preserving
the supported single LESION auxiliary state and rejecting unknown membership
categories.

## Source comparison

The committed Java excerpt harness pins `LesionToadSegmentation.java` and the
bundled descriptor order. With fixed class/lesion terms, unit shapes, direct
relations, and zero smoothing/prior, the Python result matches all 13 class
membership bits, the shared lesion-membership bit, and the denominator bit.
This is a one-voxel formula comparison with fixed terms, not a connected
iteration or full workflow equivalence result.

Focused command:

```text
./.venv/bin/pytest -q tests/test_lesion_descriptor_membership.py tests/test_lesion_objective_membership_pass.py tests/test_lesion_continuing_iteration.py
```

Result: **8 passed, 0 failed, 0 skipped in 5.38 seconds**.

## Full suite

Full suite on tested commit `11cb482c17282589ebe8777406644670da69f189`:
**1,149 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 785.28
seconds**, peak sampled RSS **15.498 GB** under the 24 GB watchdog. See the
[full-suite report](full-pytest-2026-10-04-descriptor-membership.md).

The suite used synthetic data only. Data-term-to-membership workflow
composition, iterative stopping, end-to-end Java comparison, and real-MRI
validation remain open.
