# AtlasEM focused regression — 2026-10-16

Code under test: commit `a09d8c20` (`Port AtlasEM rigid pyramid preparation`). This includes the AtlasEM rigid-prior, crop/uncrop, registration-recrop, and registration-pyramid preparation slices.

Command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_classification.py \
  tests/test_atlas_em_gaussian_memberships.py tests/test_atlas_em_iteration.py \
  tests/test_atlas_em_preprocessing.py tests/test_atlas_em_prior_normalization.py \
  tests/test_atlas_em_rigid_priors.py tests/test_atlas_em_cropping.py \
  tests/test_atlas_em_registration_pyramid.py \
  tests/test_atlas_em_statistics.py tests/test_atlas_em_outputs.py \
  tests/test_mipav_morphology.py tests/test_plus_atlas_em_dura_integration.py
```

Result: **22 passed in 5.76 s**. This covers synthetic selected-source method excerpts, Python component behavior, the AtlasEM registration-result recrop sequence, registration-pyramid preparation, and the synthetic AtlasEM-output-to-Dura handoff. Crop-bound reduction and uncropping are vectorized for brain-resolution volumes while preserving the source-comparison results. It does not validate the missing selected `ibsr8obj` bundle, the AtlasEM rigid optimizer, registration convergence, a full PLUS run with the historical pre-stages, or fresh MATLAB/MIPAV/JIST workflow equivalence. The last full-suite result remains the separate report at [`full-pytest-2026-10-16-atlasem-dura-handoff.md`](full-pytest-2026-10-16-atlasem-dura-handoff.md).

## Follow-up after two-stage registration composition

Code under test: commit `61d884c7` (`Check selected AtlasEM demons defaults`), including the composed two-stage registration and selected plugin defaults.

Command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_descriptor.py \
  tests/test_atlas_em_descriptor_initialization.py \
  tests/test_atlas_em_registration_from_descriptor.py \
  tests/test_atlas_em_inventory_cli.py tests/test_atlas_em_iteration.py \
  tests/test_atlas_em_gaussian_memberships.py tests/test_atlas_em_cli.py \
  tests/test_atlas_em_prior_normalization.py tests/test_atlas_em_cropping.py \
  tests/test_atlas_em_classification.py tests/test_atlas_em_preprocessing.py \
  tests/test_atlas_em_outputs.py tests/test_atlas_em_rigid_priors.py \
  tests/test_atlas_em_statistics.py tests/test_atlas_em_registration_pyramid.py \
  tests/test_plus_atlas_em_dura_integration.py
```

Result: **40 passed, 3 expected RAS orientation warnings in 6.11 s**. This
includes the descriptor-driven two-stage rigid-registration/recrop composition
and the shared source-faithful rigid-update Java fixture. This remains
synthetic and selected-slice evidence; it does not supply the missing
`ibsr8obj` assets or establish direct `DemonToadDeformableAtlas` optimizer
parity or complete PLUS source equivalence.

## Follow-up after AtlasEM internal-mask construction

Code under test: commit `0ca74cd6` (`Build AtlasEM internal mask from image channels`).

Command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_descriptor.py \
  tests/test_atlas_em_descriptor_initialization.py \
  tests/test_atlas_em_registration_from_descriptor.py \
  tests/test_atlas_em_inventory_cli.py tests/test_atlas_em_iteration.py \
  tests/test_atlas_em_gaussian_memberships.py tests/test_atlas_em_cli.py \
  tests/test_atlas_em_prior_normalization.py tests/test_atlas_em_cropping.py \
  tests/test_atlas_em_classification.py tests/test_atlas_em_preprocessing.py \
  tests/test_atlas_em_outputs.py tests/test_atlas_em_rigid_priors.py \
  tests/test_atlas_em_statistics.py tests/test_atlas_em_registration_pyramid.py \
  tests/test_plus_atlas_em_dura_integration.py tests/test_mipav_morphology.py
```

Result: **43 passed, 3 expected RAS orientation warnings in 6.50 s**. This adds
the AtlasEM positive-channel internal mask constructor to the two-stage
registration, EM, output, and Dura synthetic coverage. It remains source-slice
evidence, not selected-atlas or full-workflow parity.

## Follow-up after registered fixed-prior execution composition

Code under test: commit `d800db56` (`Compose AtlasEM registration with fixed-prior EM`).

Command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_descriptor.py \
  tests/test_atlas_em_descriptor_initialization.py \
  tests/test_atlas_em_registration_from_descriptor.py \
  tests/test_atlas_em_inventory_cli.py tests/test_atlas_em_iteration.py \
  tests/test_atlas_em_gaussian_memberships.py tests/test_atlas_em_cli.py \
  tests/test_atlas_em_prior_normalization.py tests/test_atlas_em_cropping.py \
  tests/test_atlas_em_classification.py tests/test_atlas_em_preprocessing.py \
  tests/test_atlas_em_outputs.py tests/test_atlas_em_rigid_priors.py \
  tests/test_atlas_em_statistics.py tests/test_atlas_em_registration_pyramid.py \
  tests/test_plus_atlas_em_dura_integration.py tests/test_mipav_morphology.py
```

Result: **43 passed, 3 expected RAS orientation warnings in 5.98 s**. The
descriptor test now composes registration, recropping, prior transformation,
the initial Atlas membership update, a one-iteration fixed-prior EM loop, and
Dura-input generation. This stays within synthetic selected-slice evidence;
the one-iteration cap does not establish convergence or Dura evolution.
