# LesionTOADS covariance synthetic reference

This asymmetric three-channel RAS NIfTI phantom pins the fixed inputs for the
selected `LesionToadSegmentation.computCovariance()` method excerpt. It
contains hard-label classes 0–4 with class 5 absent, generic OBJECT,
LESION, BLACKHOLE and OUTLIER types, membership ties and zero membership,
non-collinear channel residuals, and a spatially varying multiplicative field.
The selected method traverses the full grid, including all faces.

`java_reference.json` stores bitwise class and shared LESION covariance
outputs with and without fields. The Java harness source hash is included in
the file. Regenerate the phantom with
`python scripts/generate_lesion_covariance_fixture.py` and its selected-method
reference with `scripts/generate_java_lesion_covariance_reference.sh`.

This is a dormant statistics API comparison. The selected wrapper's call to
`computCovariance()` is commented out; this port does not enable it or claim
that covariance participates in the active workflow.
