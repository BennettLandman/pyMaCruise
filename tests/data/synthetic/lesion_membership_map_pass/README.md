# LesionTOADS membership map pass fixture

This deterministic asymmetric, folded synthetic input exercises the supported
map-update composition for `OBJECT + LESION + BLACKHOLE + OUTMASK`. The image
`inputs_ras.nii.gz` is a float32 4-D NIfTI with shape `7 x 8 x 9 x 24`, a
positive RAS affine, and channels in this order:

1. Initial class memberships 0–4.
2. Initial lesion membership, black-hole membership, and outlier membership.
3. Class data terms 0–4, then lesion, black-hole, and OUTMASK data terms.
4. Shape values 0–4, including an asymmetric folded support and two local
   zero-shape holes.
5. Ventricle, GM, and brainstem distance maps encoded as float32 integers.

The Python fixture generator records the input construction. The selected-Java
reference harness receives an x/y/z/channel text export of those exact NIfTI
float32 values, then performs the source-ordered interior membership update
and emits every interior class/auxiliary map value as IEEE float bits. The
Python test compares the complete interior outputs bitwise and checks that all
six boundary faces remain unchanged. Regenerate with:

```sh
.venv/bin/python tests/data/synthetic/lesion_membership_map_pass/generate_fixture.py
.venv/bin/python tests/data/synthetic/lesion_membership_map_pass/generate_reference.py
```

The reference is pinned to the selected VandyCRUISE
`LesionToadSegmentation.java` source SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.

This is a full-map membership composition reference for the supported class
types and supplied fixed data terms. It does not validate statistics
estimation, data-term construction, energy/objective accumulation, unsupported
object types, distance-map creation, or the complete LesionTOADS workflow.
