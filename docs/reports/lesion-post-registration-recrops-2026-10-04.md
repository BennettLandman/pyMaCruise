# LesionTOADS post-registration recrop — 2026-10-04

## Scope

P39zzzzzzbk composes the selected wrapper sequence after atlas registration:
`uncropImages` using the preprocessor transform, topology-aware
`findCroppingBoundaries`, `cropImages`, update atlas dimensions, sample atlas
topology and shape maps on the recropped grid using the separately optimized
atlas transform, and initialize constructor label state. The cropped image
arrays and atlas maps are in RAS voxel-axis order.

Normalization, estimation of the initial preprocessing transform, transformed-
template centroid refresh, segmentation memberships, and thinning are still
outside this slice. The resulting constructor state is not a whole LesionTOADS
run or a composite Java-output match.

## Synthetic Java reference

The test uses the committed asymmetric two-channel RAS NIfTI phantom in
[`lesion_preprocessing_crop_rigid`](../../tests/data/synthetic/lesion_preprocessing_crop_rigid/).
The selected Java excerpt harness performs first crop → uncrop → second crop.
It records inclusive bounds, dimensions, float32 transform bits, and SHA-256
digests of complete channel-first crop arrays. Its source and generated JSON
hashes are pinned by `tests/test_lesion_registered_recrop.py`; regeneration is
documented in the fixture README. This is selected-method Java evidence, not
execution in the full MIPAV/JIST runtime.

Focused command:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q \
  tests/test_lesion_registered_recrop.py \
  tests/test_lesion_initial_alignment_bounds.py \
  tests/test_lesion_preprocessing_crop.py \
  tests/test_lesion_preprocessing_crop_rigid.py \
  tests/test_lesion_preprocessing_uncrop.py
```

Result: **6 passed, 0 failed, 0 skipped in 7.03s** on implementation commit
`d3096ea6c4cb403af1dd6034f36ae002c6c576ea`.

The test confirms crop bounds/dimensions, preprocessing-transform bits, and
every crop voxel against the Java excerpt reference. It also uses a different
atlas transform and proves that it affects transformed atlas data but does not
seed preprocessing crop bounds. Atlas shape maps and constructor labels are
checked on the recropped grid; input arrays are unchanged. A left-to-right
topology input emits the RAS warning and is reoriented before calculation.
The resampling and constructor components are compared with direct Python
component calls; no composite Java segmentation comparison is claimed.

## Full-suite evidence

Command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Tested implementation commit: **`d3096ea6c4cb403af1dd6034f36ae002c6c576ea`**.
Result: **1,161 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
880.37s**. Peak sampled RSS was **12.170 GB**; the watchdog did not terminate
the run.

This validates synthetic regression only. Whole LesionTOADS iteration,
MaCRUISEplus, source-workflow acceptance, and human MRI validation remain open.
