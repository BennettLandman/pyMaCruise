# Registered TOADS startup to correction-field fit (P18cyav)

## Ported behavior

`fit_registered_toads_correction_fields` consumes one
`ToadsAtlasAlignmentRun`, extracts the recropped image, exact initial
membership maps, constructor segmentation, and class centroids, initializes a
unity field per channel, and invokes the selected exact global Chebyshev fit
in source channel order. The selected default profile starts at degree 3;
degree remains an explicit supported option for degree 1–3 method checks.

This is the wrapper handoff after registration and segmentation initialization.
It does not apply the fields to image intensities, rerun membership/relations,
or continue a TOADS iteration.

## Committed synthetic inputs and Java comparison

The test uses committed RAS `.nii.gz` files under
[`tests/data/synthetic/toads_registered_field_refit/`](../../tests/data/synthetic/toads_registered_field_refit/):
an asymmetric anisotropic T1 and brain mask, an 11-class atlas topology, and
one shape prior per class. The test executes selected atlas registration and
crop/constructor initialization before fitting. It then supplies the
run-derived images, memberships, labels, and centroids to the Java formula
harness for the exact degree-three field method.

The normal matrix and right-hand side match the Java formula harness exactly.
The complete output field is within `2e-6` absolute error, is exactly one
outside the correction mask, and follows the selected one-channel order.
This is synthetic integration around a selected method reference, not a
complete Java wrapper execution or full TOADS validation.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_field_refit.py tests/test_toads_registered_atlas_membership_integration.py tests/test_toads_exact_polynomial_field_degrees23.py
```

Result: **13 passed, 0 failed, 0 skipped in 30.21 seconds** on implementation
and fixture commit `a4c0a581de32eb3922f67cd0c9e67a607459de24`.

The full repository suite has not yet been rerun on this commit. P18cyav
remains Started pending the next major integration checkpoint.
