# Full pytest evidence: P39zzo constrained LesionTOADS centroids

- Run: `.venv/bin/pytest -q`
- Result: **609 passed, 2 warnings, in 79.97s**.
- Warnings: the two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `68dab53f1f851dcc46053fc4c0c317c6cff28f07`.
- Focused selected-Java comparison: `tests/test_lesion_constrained_centroids.py` (**2 passed**).
- Synthetic inputs and reference provenance: [constrained-centroid phantom](../../tests/data/synthetic/lesion_constrained_centroids/README.md).

The Python function and selected-method Java harness match raw float32 output
bits for both uncorrected and field-corrected centroid updates over two
channels and four classes. The case exercises the source's 45% upper bounds,
the empty-class zero mean followed by its lower bound, and its `-1e10` first
lower-bound sentinel. Accumulation follows source channel/class and x/y/z
order. Invalid grid, membership-class, centroid-shape, and field-shape inputs
are also checked.

This validates `computeConstrainedCentroids()` on fixed synthetic labels and
memberships; it does not validate the broader class/lesion/black-hole
centroid, variance, channel-weight, membership, or LesionTOADS iteration
pipeline.
