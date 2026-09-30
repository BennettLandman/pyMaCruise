# Full pytest evidence: P39zzs LesionTOADS pooled noise variance

- Run: `.venv/bin/pytest -q`
- Result: **616 passed, 2 warnings, in 113.03s**.
- Warnings: the two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `710085d2fe9407a27c77fcd641c17fa476c5ac33`.
- Focused selected-Java pooled-variance tests: `tests/test_lesion_noise_variance.py` (**2 passed**).
- Synthetic input and reference provenance: [LesionTOADS variance phantom](../../tests/data/synthetic/lesion_variances/README.md).

Python and the selected-method Java excerpt match raw float32 standard-
deviation bits on the phantom with and without field correction. The update
pools membership-squared residuals over hard classes 1 through N−1, excludes
class 0, and writes the same per-channel standard deviation to every class.
The suite checks class-0 intensity invariance, an empty non-background pool
that returns zero for every class, and malformed centroid input.

This validates the selected pooled noise estimator on fixed inputs only; it
does not validate the enclosing statistics or membership-update iteration.
