# AtlasEM rigid-registration and recrop handoff — 2026-10-16

Added `register_and_recrop_atlas_em` to connect the new AtlasEM rigid-pyramid adapter with the selected preprocessing recrop sequence. It passes the registration's final transform into the source-checked sequence: install the crop-space transform, uncrop channels and transform to native coordinates, recompute the foreground bounds, then crop with the configured zero border.

The regression uses deterministic synthetic 12³ data, a crop with a two-voxel halo, three membership/prior classes, one rigid-registration iteration, and two pyramid levels. It checks that the registration produces a finite transform and that recropping returns finite channel data on the original native grid. This is a functional Python synthetic handoff; it does not validate the selected Java optimizer's final transform or the source-configured AtlasEM initialization.

Focused command on implementation commit `99de736c0900519a8bee27c61fdcf3a894c12ffa`:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_classification.py tests/test_atlas_em_gaussian_memberships.py \
  tests/test_atlas_em_iteration.py tests/test_atlas_em_preprocessing.py \
  tests/test_atlas_em_prior_normalization.py tests/test_atlas_em_rigid_priors.py \
  tests/test_atlas_em_cropping.py tests/test_atlas_em_registration_pyramid.py \
  tests/test_atlas_em_statistics.py tests/test_atlas_em_outputs.py tests/test_mipav_morphology.py \
  tests/test_plus_atlas_em_dura_integration.py
```

Result: **24 passed, 0 failed, 0 skipped in 5.73 s**. The selected Java pyramid-preparation excerpt remains the only Java-backed evidence in this focused command. Full AtlasEM registration parity, selected atlas assets, and complete PLUS integration remain open. No full-suite run was performed; the TOADS source-method line estimate is unchanged.
