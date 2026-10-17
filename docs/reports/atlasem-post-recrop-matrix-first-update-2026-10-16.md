# AtlasEM post-recrop matrix-first rigid update — 2026-10-16

Implementation commit: `8685672cc6121854185bbd3e1b3d63acddb2b9bc`

## Source finding

The selected `DemonToadDeformableAtlas` class sets `precompute=true` for a
linear transform in its initialization path and builds the image-to-template
matrix before coefficient accumulation. Its exact-registration loop then
uses the matrix-first coordinate helper; the post-recrop stage retains that
configuration. The class source SHA-256 is
`96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a`.

The AtlasEM post-recrop adapter now selects `precomputed_affine=True` in the
shared Python rigid kernel. Before this change it followed the direct
coordinate path. This aligns the float32 coordinate evaluation route with the
selected linear Java configuration. It does not prove direct Java optimizer
output parity: registration coefficients, trial path and resulting transform
have not yet been compared against a run of the complete
`DemonToadDeformableAtlas` class.

## Validation

The adapter test compares the complete one-iteration transform, energy,
lambda, relative change, accepted/rejected decision and trial count against
the shared kernel with the selected offset-zero, stride-two, source-faithful,
matrix-first settings. The AtlasEM/registration/Dura/morphology group, with
the Toads rigid-update reference test included, passed **47 tests, 3 expected
RAS-orientation warnings, in 7.01 seconds** on the implementation commit.

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_descriptor.py \
  tests/test_atlas_em_descriptor_initialization.py \
  tests/test_atlas_em_registration_from_descriptor.py \
  tests/test_atlas_em_inventory_cli.py tests/test_atlas_em_iteration.py \
  tests/test_atlas_em_gaussian_memberships.py tests/test_atlas_em_cli.py \
  tests/test_atlas_em_cli_workflow.py tests/test_atlas_em_prior_normalization.py \
  tests/test_atlas_em_cropping.py tests/test_atlas_em_classification.py \
  tests/test_atlas_em_preprocessing.py tests/test_atlas_em_outputs.py \
  tests/test_atlas_em_rigid_priors.py tests/test_atlas_em_statistics.py \
  tests/test_atlas_em_registration_pyramid.py \
  tests/test_plus_atlas_em_dura_integration.py tests/test_mipav_morphology.py \
  tests/test_atlas_em_synthetic_full_handoff.py \
  tests/test_toads_rigid_registration_update.py
```

This is selected-source configuration evidence plus a shared-kernel
comparison. It is not a compiled execution of the selected Java class or
selected-atlas behavior. The unavailable `ibsr8obj` payload and direct Java
optimizer comparison remain open.
