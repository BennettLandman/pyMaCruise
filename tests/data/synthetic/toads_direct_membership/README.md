# TOADS direct membership-update phantom

This deterministic synthetic fixture contains three 7×6×5×4 float32 RAS
`.nii.gz` arrays on a 1 mm isotropic affine: per-class input data energy,
objective energy, and initial memberships. Its generator defines every value
and inserts energies at zero, exactly at the source `ZERO` threshold, just
below it, and just above it.

`java_reference.json` records direct optimized-class reciprocal-energy
conversion and class normalization from a selected
`MultipleToadSegmentation.computeExactMemberships` formula excerpt. It also
records per-class objective contributions after normalization. The source
float32 constants are `ZERO=1e-10f` and `INF=1e10f`; the selected source file
SHA-256 is included in the JSON. The harness intentionally exercises direct
relations with smoothing and shape-prior contributions set to zero.

Regenerate the NIfTIs with `scripts/generate_toads_direct_membership_fixture.py`
and the Java formula reference with
`scripts/generate_java_toads_direct_membership_reference.sh`.
