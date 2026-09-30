# Full pytest evidence: P39zzv optimized LesionTOADS centroid weights

- Run: `.venv/bin/python -m pytest -q`
- Result: **622 passed, 2 warnings, in 113.80s**.
- Warnings: two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `92d6eb71e420986a404442c2b7a4dcdb6e3bce76`.
- Focused selected-Java tests: `tests/test_lesion_optimized_centroids.py` (**2 passed**).
- Synthetic input and reference provenance: [optimized-centroid phantom](../../tests/data/synthetic/lesion_optimized_centroids/README.md).
- PowerTable provenance: preserved selected source `PowerTable.java`, SHA-256 `2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae`.

On the 5×4×3 RAS phantom, raw float32 centroid bits match for `none`,
`linked`, `damped`, `prior`, and fallback updates, with and without fields.
The selected-method harness invokes the preserved Java PowerTable. Cases
cover optimized CSF/GM/WM classes, CSF's T1_MPRAGE versus FLAIR fallback,
per-channel exponents, an exact zero residual, a residual below the `ZERO`
cutoff, cached table values, direct lookup outside the cache range, and an
absent class.

This validates optimized tissue-centroid weights for supplied maps only; it
does not validate optimized LESION/BLACKHOLE branches or full iteration.
