# P18al: optimized-tissue damped centroid update

## Scope

This slice implements the selected `MultipleToadSegmentation.computeExactCentroids`
update for the optimized CSF, GM, and WM tissue classes. It computes the
membership-squared weighted centroid update using the source `PowerTable`
weighting and then applies the wrapper-selected damped centroid mode with
`centroidSmoothness=0.1f`. The optional intensity field multiplies image
intensities before residual calculation.

This is one update primitive. It is not a complete TOADS segmentation, does
not validate convergence or registration, and does not cover object, lesion,
outlier, Rayleigh, linked-centroid, or prior-centroid modes.

## Reference and fixtures

The committed fixture reuses RAS, 1 mm compressed-NIfTI images and a field
from `tests/data/synthetic/toads_optimized_likelihood/`, along with the first
three membership channels from `tests/data/synthetic/toads_direct_membership/`.
The shape is 7×6×5 with three modalities and three tissue classes. The test
uses initial centroids and exponents stored in
`tests/data/synthetic/toads_optimized_centroid/java_reference.json`.

The Java harness calls the selected `PowerTable.java` class directly. Its
recorded SHA-256 is
`2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae`; the
selected `MultipleToadSegmentation.java` SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
The harness source/formula and exact output arrays can be regenerated with
`scripts/generate_java_toads_optimized_centroid_reference.sh`.

## Validation

`tests/test_toads_optimized_centroids.py` compares all nine CSF/GM/WM modality
centroids against Java for both field-disabled and field-enabled modes using
exact float32 array equality. The fixture includes an image sample exactly
equal to an initial centroid, exercising the source's zero-distance
`PowerTable` factor fallback. It also checks empty class support leaves that
centroid unchanged and rejects invalid array shapes. Six focused tests pass.

The complete configured suite passed **323 tests in 50.84 seconds** on
implementation commit
`3b2f7c12a98ad718f4de1f98fe37ab2b894f47e9`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md). This matches the Java
reference for the tested formula slice and does not establish full TOADS
equivalence.
