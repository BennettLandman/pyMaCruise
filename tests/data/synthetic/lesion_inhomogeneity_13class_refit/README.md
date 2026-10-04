# Synthetic 13-class LesionTOADS field-refit fixture

This deterministic `(21, 21, 21)` RAS phantom starts from the committed
`toads_exact_polynomial_field_degrees23` synthetic image, memberships,
segmentation, field, and centroids. It adds the two LesionTOADS descriptor
classes omitted from that base profile: `lesion` and a thirteenth `obj`.
Added class locations and memberships are deterministic. The complete class
type sequence is in `class_types.json` and follows the bundled 13-class
LesionTOADS descriptor order.

All images are compressed NIfTI and retain the source phantom's RAS affine.
Regenerate from the repository root with:

    ./.venv/bin/python tests/data/synthetic/lesion_inhomogeneity_13class_refit/generate_fixture.py

The fixture validates the selected exact polynomial field kernel for a
13-class layout and the two-channel fitting adapter against the pinned Java
formula harness. It is synthetic component evidence, not full LesionTOADS or
MaCRUISEplus workflow validation.
