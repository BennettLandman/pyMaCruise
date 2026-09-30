# Full pytest evidence: P39zzw LESION/BLACKHOLE centroid routing

- Run: `.venv/bin/python -m pytest -q`
- Result: **624 passed, 2 warnings, in 155.83s**.
- Warnings: two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `62d269854156c7fc45b2b35ab8f6c851e95b4e29`.
- Focused selected-Java tests: `tests/test_lesion_auxiliary_centroids.py` (**2 passed**).
- Synthetic input and reference provenance: [LESION/BLACKHOLE centroid phantom](../../tests/data/synthetic/lesion_special_centroids/README.md).
- PowerTable provenance: preserved selected source `PowerTable.java`, SHA-256 `2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae`.

The selected source excerpt and Python port match raw float32 class, lesion,
and black-hole centroids in all five class-centroid modes, with and without
field correction. The phantom routes 14 hard-labeled samples to class
accumulators, 7 to lesion accumulators, and 9 to black-hole accumulators. It
includes class/lesion, lesion/black-hole, and all-three ties; the exact source
tie precedence is preserved. Lesion updates use PowerTable weights and the
selected 20% blend, while black-hole updates use the direct selected mean.
An absent class and malformed map/parameter cases are covered.

This validates centroid routing for fixed membership maps only; it does not
estimate auxiliary memberships or validate full LesionTOADS iteration.
