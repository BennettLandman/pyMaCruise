# LesionTOADS descriptor membership entry (started)

The raw-energy adapter now has a downstream seam for the bundled descriptor.
This slice maps its ordered topology labels to the class categories used by
the existing in-place membership updater: `obj` to OBJECT, `mask` to MASK,
`optimized_csf` to the ordinary MASK-like membership branch, and `lesion` to
LESION. Optimized-CSF's raw term remains the source's active outlier fallback,
implemented upstream in P39zzzzzzbd.

The generic updater previously required one BLACKHOLE and one OUTMASK class.
Those classes are absent from the bundled 13-class LesionTOADS descriptor, so
validation now allows those auxiliary families to be omitted while retaining
the existing exact-one LESION requirement for descriptor composition. Other
class families remain explicitly rejected by this descriptor adapter.

The committed Java excerpt harness now pins the selected source SHA and
descriptor class order. With fixed raw terms, unit shapes, direct relations,
and zero smoothing/prior, all 13 class-membership bits, the shared lesion
membership bit, and denominator bit match the Java harness. This is a bounded
single-voxel formula comparison; it does not validate data-term-to-membership
workflow composition or a connected iteration. Command:

```text
./.venv/bin/pytest -q tests/test_lesion_descriptor_membership.py tests/test_lesion_objective_membership_pass.py tests/test_lesion_continuing_iteration.py
```

Result: **8 passed, 0 failed, 0 skipped in 5.38 seconds**. The 3×3×3
deterministic case exercises all 13 descriptor classes at one interior voxel,
checks class plus lesion normalization, and confirms boundary state is left
untouched. Full-suite validation remains pending before marking the plan row
complete.
