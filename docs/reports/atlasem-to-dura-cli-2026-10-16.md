# AtlasEM to pre-LesionTOADS Dura command — 2026-10-16

Implementation commit: `b15e0070b890ac03867209dff87f28e73d86c4e0`

## Runnable interface

```sh
pymacruise atlasem-run \
  --t1 subject_t1.nii.gz \
  --brain-mask subject_brain_mask.nii.gz \
  --descriptor atlas/descriptor.txt \
  --output-dir output/atlasem-plus \
  --provenance output/atlasem-plus/provenance.json
```

The command requires a binary brain mask on the same physical RAS grid as the
T1. It runs the descriptor-driven single-channel AtlasEM path through both
registration stages, fixed-prior EM and native-grid output persistence, then
passes the saved memberships and supplied mask into the selected
pre-LesionTOADS Dura adapter. The Dura iteration cap defaults to five and can
be lowered for synthetic runs with `--max-dura-iterations`.

The command writes 13 compressed RAS NIfTIs: memberships, tissue-index
segmentation, AtlasEM's internal mask, AtlasEM Dura-input channels,
pre-LesionTOADS Dura-input channels and segmentation, stripped T1, corrected
white/gray maps, brain-mask level set and mask, and WM/GM topology products.
Optional provenance records input paths, settings, EM/Dura completed iteration
counts, and output paths.

## Validation

The synthetic CLI fixture creates a 12³ T1, matching binary brain mask,
three-class descriptor, and all compressed NIfTI atlas assets. It runs one
registration iteration per stage, one EM iteration, and one Dura pass, then
checks the complete output set for finite voxels and RAS geometry and confirms
the provenance iteration counts.

Focused AtlasEM/registration/Dura/morphology regression: **45 passed, 3
expected orientation warnings, in 6.57 seconds** on implementation commit
`b15e0070b890ac03867209dff87f28e73d86c4e0`. The runnable help was inspected
through the installed CLI entry point.

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

This is synthetic workflow execution, not selected-atlas validation, direct
`DemonToadDeformableAtlas` optimizer parity, natural human convergence, or
whole-PLUS source equivalence. The command currently has optional provenance
but does not yet support stage-level hash-verified resume. Obtaining the
selected `ibsr8obj` descriptor and all referenced atlas files remains needed
for source-configured PLUS execution.
