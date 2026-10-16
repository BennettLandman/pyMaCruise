# AtlasEM rigid-registration adapter — 2026-10-16

Added `register_atlas_em_rigid_pyramid` to connect the AtlasEM rigid-registration profile to the existing source-compared transform sampling, coefficient, line-search, and pyramid routines.

The adapter follows the selected `AlgorithmAtlasEMSegmentation` and `DemonToadDeformableAtlas` configuration: stride 2, zero offset, lambda 1 at each scale, lambda factor 1.5, rejection limit 10, relative-change cutoff `1e-6`, minimum lambda `0.001`, and no scale-1 optimization. It selects `source_faithful=False` because `DemonToadDeformableAtlas.computeScaledRegistrationCoefficients` uses the direct gradient divided by the prior norm. The adjacent `MultipleToadDeformableAtlas` has a distinct historical Hessian/gradient mix-up and is not the selected AtlasEM class.

The adapter uses the previously validated TOADS Python kernels and deliberately does not imply whole-class equivalence. The selected AtlasEM Java optimizer has not yet been compared end-to-end on the same synthetic arrays; its 8-object descriptor and atlas payload are also unavailable. AtlasEM registration parity and asset-backed initialization remain open.

Focused command on implementation commit `7983f03ae15ff025a74c4a6e4c2939de829ebcb3`:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_registration_pyramid.py tests/test_atlas_em_rigid_priors.py tests/test_atlas_em_cropping.py
```

Result: **4 passed, 0 failed, 0 skipped in 2.03 s**. The selected Java-excerpt check in this run validates pyramid preparation; the new adapter test verifies the source-profile dispatch and result handoff against the configured shared Python kernel. This is synthetic adapter evidence, not a direct Java `DemonToadDeformableAtlas` comparison. No full-suite run was performed, and no TOADS source-line estimate changed.
