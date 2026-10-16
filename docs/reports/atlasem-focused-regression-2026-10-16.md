# AtlasEM focused regression — 2026-10-16

Code under test: commit `341c8412` (`Connect AtlasEM registration recropping stages`). This includes the AtlasEM rigid-prior, crop/uncrop, and registration-recrop slices.

Command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_classification.py \
  tests/test_atlas_em_gaussian_memberships.py tests/test_atlas_em_iteration.py \
  tests/test_atlas_em_preprocessing.py tests/test_atlas_em_prior_normalization.py \
  tests/test_atlas_em_rigid_priors.py tests/test_atlas_em_cropping.py \
  tests/test_atlas_em_statistics.py tests/test_atlas_em_outputs.py \
  tests/test_mipav_morphology.py tests/test_plus_atlas_em_dura_integration.py
```

Result: **21 passed in 5.28 s**. This covers synthetic selected-source method excerpts, Python component behavior, the AtlasEM registration-result recrop sequence, and the synthetic AtlasEM-output-to-Dura handoff. Crop-bound reduction and uncropping are vectorized for brain-resolution volumes while preserving the source-comparison results. It does not validate the missing selected `ibsr8obj` bundle, configured AtlasEM registration, a full PLUS run with the historical pre-stages, or fresh MATLAB/MIPAV/JIST workflow equivalence. The last full-suite result remains the separate report at [`full-pytest-2026-10-16-atlasem-dura-handoff.md`](full-pytest-2026-10-16-atlasem-dura-handoff.md).
