# AtlasEM rigid-registration pyramid preparation — 2026-10-16

`build_atlas_em_rigid_registration_pyramid` reuses the existing source-checked `ImageFunctions.subsample` port through `build_toads_registration_pyramid`. The selected `DemonToadDeformableAtlas.registerShapesPyramid` starts at scale `2**(levels-1)`, halves the scale to two, and does not optimize a scale-one level. It constructs both the membership volume and every class shape volume from the original arrays at each scale.

The Java excerpt calls the pinned `ImageFunctions.subsample` implementation for scale factors 4 and 2 on deterministic 8×8×8, three-class arrays. Every sampled membership and shape value matches Python bit-for-bit. The reference records `DemonToadDeformableAtlas.java` SHA-256 `96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a` and the copied Java `ImageFunctions.java` SHA-256 `be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892`.

Focused validation on code commit `a09d8c20`:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_registration_pyramid.py
```

Result: **1 passed in 1.03 s**. This validates pyramid preparation, not AtlasEM's optimizer, registration convergence, selected atlas assets, or full PLUS integration.
