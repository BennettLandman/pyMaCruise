# Full pytest evidence: P39zzu OUTLIER centroid suppression

- Run: `.venv/bin/python -m pytest -q`
- Result: **620 passed, 2 warnings, in 110.09s**.
- Warnings: two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `c87cbce7ad2023529cd64a4c018047815a34a05f`.
- Focused selected-Java tests: `tests/test_lesion_outlier_centroids.py` (**2 passed**).
- Synthetic input and reference provenance: [OUTLIER centroid phantom](../../tests/data/synthetic/lesion_outlier_centroids/README.md).

For classes declared OUTLIER, both the selected Java excerpt and Python port
replace membership by `membership × (1 − outlierMembership)` before squaring
and accumulating. Raw float32 centroid bits match for `none`, `linked`,
`damped`, `prior`, and fallback modes, with and without field correction. The
map covers 0, 0.25, 0.5, 0.75, and 1; controls confirm OBJECT/OUTMASK types
remain generic. Missing maps, wrong map shapes, and nonfinite values reject.

This validates the OUTLIER centroid accumulation for fixed maps only; it does
not estimate outlier membership or validate the complete LesionTOADS loop.
