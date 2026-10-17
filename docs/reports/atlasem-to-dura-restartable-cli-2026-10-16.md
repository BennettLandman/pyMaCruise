# Restartable AtlasEM-to-Dura synthetic workflow — 2026-10-16

Implementation commit: `a9ced406f4f802980bef0e26ebb4c3d1e5964e2e`

## Runnable command

```sh
pymacruise atlasem-run \
  --t1 subject_t1.nii.gz \
  --brain-mask subject_brain_mask.nii.gz \
  --descriptor atlas/descriptor.txt \
  --output-dir output/atlasem-plus

# Continue a completed or interrupted run after verifying its inputs and outputs.
pymacruise atlasem-run \
  --t1 subject_t1.nii.gz \
  --brain-mask subject_brain_mask.nii.gz \
  --descriptor atlas/descriptor.txt \
  --output-dir output/atlasem-plus \
  --resume
```

Every run records `output/atlasem-plus/run.json`. It fingerprints the T1, the
T1-grid binary brain mask, descriptor and every referenced atlas asset, the
effective settings, and the Python implementation. It records independent
AtlasEM and pre-LesionTOADS Dura stages. Each stage record includes output
paths, SHA-256, byte size, dimensions and RAS axis codes. Resume requires an
identical run signature and verifies the expected output paths and hashes.

When AtlasEM is intact and a Dura output is damaged, resume reloads the
verified native-grid memberships and reruns Dura only. If a registration or
EM output is damaged, AtlasEM reruns and the Dura stage is recomputed against
the new AtlasEM output signature. A changed input, atlas, setting, or Python
implementation rejects resume. Use `--overwrite` for a fresh run in an
existing output directory.

## Synthetic evidence

The 12³ test creates a three-class atlas descriptor and compressed payloads,
then runs one iteration at each registration stage, one EM iteration and one
Dura pass. It checks all 13 product NIfTIs for finite data and native RAS
geometry, and verifies the manifest stage order and iteration counts. It then
resumes unchanged and confirms that history does not grow; corrupts one Dura
image and confirms the AtlasEM history remains unchanged while Dura reruns;
and edits the brain mask to confirm signature mismatch prevents reuse.

The AtlasEM/registration/Dura/morphology focused regression passed **45
tests, 3 expected RAS orientation warnings, in 6.98 seconds** on the
implementation commit above.

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

This demonstrates synthetic command execution, stage validation and resume.
It is not selected `ibsr8obj` execution, direct
`DemonToadDeformableAtlas` optimizer parity, natural human convergence, or
whole-PLUS source equivalence. The selected `ibsr8obj` payload remains
unavailable. The AtlasEM and Dura iteration limits in this phantom are
deliberately small and do not represent a clinical/default runtime estimate.
