# Descriptor-driven AtlasEM command — 2026-10-16

Implementation commit: `c5b1856489d73ae1b5f654fb026c2ccb1a9aa67a`

## Runnable interface

```sh
pymacruise atlasem-run \
  --t1 subject_t1.nii.gz \
  --descriptor atlas/descriptor.txt \
  --output-dir output/atlasem \
  --provenance output/atlasem/run.json
```

The command accepts one 3D T1 image and a descriptor whose referenced
registration assets are `.nii.gz`. It runs robust intensity normalization,
source-shaped crop preparation, descriptor intensity initialization, the two
registered rigid stages with the intervening native-grid recrop, fixed-prior
EM, and native-grid RAS restoration. The selected plugin defaults are used
unless iteration controls are overridden. It writes memberships, the
tissue-index segmentation, AtlasEM's internal mask, and the three-channel
Dura-input image. The optional provenance JSON records parameters and paths.

The command uses the supplied AtlasEM internal-mask rule (positive image
channels followed by MIPAV hole filling). It does not consume or replace the
outer MaCRUISE brain mask. Dura can be run after this command with the existing
`atlasem-prelesion-dura` interface using the saved memberships and the user's
T1-grid brain mask.

## Validation

The synthetic test constructs a 12³ T1, a three-class descriptor, and all four
compressed NIfTI atlas assets in a temporary directory, then invokes the CLI dispatcher with bounded iteration counts; the installed
entry point help was also inspected. It verifies all four
outputs, native RAS orientation, finite voxels, and the provenance summary.
The broader AtlasEM/registration/Dura/morphology slice passed **45 tests, with
3 expected RAS-orientation warnings, in 6.49 seconds** on the implementation
commit above. The CLI `--help` output was also inspected.

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
  tests/test_atlas_em_synthetic_full_handoff.py
```

This validates command wiring and synthetic execution. It is not evidence of
selected `ibsr8obj` behavior, direct `DemonToadDeformableAtlas` optimizer
parity, natural convergence on a human T1, or complete PLUS equivalence. The
selected `ibsr8obj` descriptor and payloads remain unavailable. The command
does not yet persist intermediate restart checkpoints or execute Dura in the
same invocation.
