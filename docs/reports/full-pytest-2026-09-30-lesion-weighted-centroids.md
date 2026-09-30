# Full pytest evidence: P39zzp weighted LesionTOADS centroids

- Run: `.venv/bin/pytest -q`
- Result: **610 passed, 2 warnings, in 79.91s**.
- Warnings: the two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `8cb8ac17c57a54c4ece1e5400dcfb61500f41dce`.
- Focused selected-Java centroid tests: `tests/test_lesion_constrained_centroids.py` (**3 passed**, including the constrained-centroid cases).
- Synthetic inputs and references: [centroid phantom](../../tests/data/synthetic/lesion_constrained_centroids/README.md).

The Python implementation and selected-method Java reference match raw
float32 output bits for factors 0, 0.35, and 1, both with and without field
correction. The test checks that factor 1 retains centroids for populated
classes, that an empty class is set to zero regardless of factor, and that
input arrays are not modified. The update uses membership-squared weights and
source channel/class/x/y/z accumulation order.

This validates `computeWeightedCentroids(factor)` on fixed hard labels and
memberships. It does not compose the broader lesion/black-hole centroid,
variance, channel-weight, or membership-estimation stages.
