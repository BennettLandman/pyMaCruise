# Full pytest evidence: P39zzx optimized LESION/BLACKHOLE centroids

- Run: `.venv/bin/python -m pytest -q`
- Result: **626 passed, 2 warnings, in 111.20s**.
- Warnings: two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `9c785fa309bb38889736d6b494e0273e50182d87`.
- Focused optimized selected-Java tests: `tests/test_lesion_optimized_auxiliary_centroids.py` (**2 passed**); ordinary auxiliary-centroid regression tests also pass (**2 passed**).
- Synthetic input and reference provenance: [optimized auxiliary-centroid phantom](../../tests/data/synthetic/lesion_optimized_auxiliary_centroids/README.md).
- PowerTable provenance: preserved selected source `PowerTable.java`, SHA-256 `2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae`.

The Java excerpt and Python port match raw float32 class, lesion, and
black-hole centroids for `OPTIMIZED_LESION` and `OPTIMIZED_BLACKHOLE` in all
five update modes, with and without fields. The selection counts are 16 class,
7 lesion, and 7 black-hole routes, with class/lesion, lesion/black-hole, and
all-three ties. The factor is computed from the winning membership squared
times corrected intensity, then applied to the selected membership-squared
accumulator. The fixture covers a factor below `ZERO=1e-10` and direct
PowerTable evaluation outside its cached maximum.

This validates optimized auxiliary centroid routing for fixed maps only; it
does not estimate lesion/black-hole maps or validate full iteration.
