# P18cyaf: signed RAS spacing handoff

Implementation commit: `0d15a143d4f8321266e82e8272c5910bac92fb56`.

The startup-to-registration adapter previously passed the absolute value of
`ToadsRASGeometry.signed_image_resolutions_mm`. Selected
`MultipleToadPreprocess.getSignedImageResolutions()` carries RAS axis direction
as `(-sx, -sy, +sz)`, and `ParametricTransform.imageToTemplate` uses those signed
values for physical coordinates. The adapter now preserves that tuple in the
inspectable setup and passes it to the rigid pyramid. Descriptor template
spacing remains independently supplied.

Source pin: `MultipleToadPreprocess.java` SHA-256
`31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7`.
Registration mapping source pins and landmarks are in
`tests/test_toads_registration_sampling.py` and its committed Java reference.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_alignment_startup.py tests/test_toads_registration_sampling.py tests/test_toads_registration.py
```

Result on the exact implementation commit: **43 passed, 0 failed, 0 skipped in
1.25 seconds**. This includes `test_ras_signed_spacing_landmarks_match_selected_java_mapping` and
an adapter assertion for `(-1.0, -1.0, 1.0)` on a unit RAS grid.

No subject MRI or full-suite run was used. The change repairs a source-coordinate
input contract; it does not validate the complete Java registration wrapper or
TOADS workflow.
