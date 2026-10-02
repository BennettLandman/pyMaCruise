# P18cyae: selected TOADS registration setup

## Source mapping

The setup record follows the selected source sequence in
`AlgorithmMultipleToads.java`: `setImageInfo`, `initShapeRegistration`,
preprocessing transform handoff, same-mode `updateShapeRegistration`, and
`registerShapesPyramid`. `initShapeRegistration` initializes a six-value rigid
transform to zero; `setTransform` assigns the preprocessing transform;
`updateShapeRegistration` retains the same memberships, sets iteration and
level counts, uses subsample 2, and rigid-to-rigid conversion returns the same
transform reference. The selected wrapper constructor sets `approx=classes`
and `useApprox=false`.

SHA-256 source pins:

- `AlgorithmMultipleToads.java`: `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`
- `MultipleToadDeformableAtlas.java`: `e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`
- `ParametricTransform.java`: `22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac`

Implementation commit: `845eecfa43b9862ce81e7d9322aeedccf4593b7b`.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_alignment_startup.py tests/test_toads_registration_sampling.py
```

Result on the exact implementation commit: **27 passed, 0 failed, 0 skipped in
1.48 seconds**. Tests verify the stored call order, membership and transform
reference identity, zero initial transform, six rigid parameters, null best map,
selected exact-mode dispatch, caller-provided limits, subsample factor, image
center, and a finite registration result from the synthetic atlas case.

## Limits

The Python setup record translates inspected source assignments; this validation
does not execute the Java wrapper or compare every internal transform-model
field against Java. It is not complete registration or TOADS workflow
equivalence. No real MRI was used, and no full-suite run was required for this
localized change.
