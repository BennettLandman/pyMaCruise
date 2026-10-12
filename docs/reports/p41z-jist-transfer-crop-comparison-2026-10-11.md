# P41z JIST transfer-plus-crop comparison — 2026-10-11

## Result

This comparison evaluates the selected legacy transfer behavior and its
immediate `CropBasedOnGM` consumer on the P41z synthetic case. The source
Dura maps are 71×89×67; T1, MaSeg, and `Atlas_GM` are 81×101×77. The source
transfer audit shows that JIST `ReorientVolume.solve(UNCHANGED,
NEAREST_NEIGHBOR)` preserves the RAS Dura arrays and their cropped extents,
then the wrapper assigns the reference header. The actual selected Java
`CubicVolumeCropper` was then run with full-grid `Atlas_GM` bounds on those
smaller maps. It returned the expected crop shape 73×91×63, padding the target
where crop indices exceeded the smaller image.

Compared with the Python route that first places Dura maps by their physical
affines on the full T1 grid and then applies the same crop, the Java crop
results differ substantially:

| Map | Differing voxels | Voxels | Maximum absolute difference | Mean absolute difference |
|---|---:|---:|---:|---:|
| CSF | 295,556 | 418,509 | 0.82494 | 0.10278 |
| GM | 282,146 | 418,509 | 0.99602 | 0.15103 |
| WM fill | 252,296 | 418,509 | 0.99628 | 0.08911 |
| WM mask | 42,457 | 418,509 | 1.0 | 0.10145 |

This establishes a source/Python discrepancy before plus surface evolution on
this cropped synthetic case. It does **not** establish that the Java result is
anatomically better, nor that this transfer behavior explains the P41z ribbon
failure. It is a selected transfer-and-crop comparison, not complete
MaCRUISEplus equivalence or real-MRI validation. The Python physical-affine
route retains the crop's physical location in the native grid and remains the
recommended behavior; source-compatible handling is an open compatibility
choice that needs an end-to-end checkpoint comparison before becoming a
workflow option.

## Reproduction

The audit combines two verified selected-source steps: the JIST reorientation
geometry audit (Java output arrays exactly equal the four source Dura arrays
in RAS order) and the direct Java `CubicVolumeCropper` replay on those arrays.
It compares those cropped outputs with Python physical transfer followed by
the same MaSeg-GM crop:

```sh
PYTHONPATH=src:scripts PYTHONNOUSERSITE=1 \
  .venv/bin/python scripts/audit_jist_reorient_native_transfer.py \
  --reference tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --map csf=build/p41z-case-normalization/iterations-0008-to-0011-surfaces/lesion-dura/lesion_toads_csf.nii.gz \
  --map gm=build/p41z-case-normalization/iterations-0008-to-0011-surfaces/lesion-dura/lesion_toads_gm.nii.gz \
  --map wm_fill=build/p41z-case-normalization/iterations-0008-to-0011-surfaces/lesion-dura/lesion_toads_wm_fill.nii.gz \
  --map wm_mask=build/p41z-case-normalization/iterations-0008-to-0011-surfaces/lesion-dura/lesion_toads_wm_mask.nii.gz \
  --output-dir build/p41z-case-normalization/iteration-0011-native-transfer-geometry-audit

PYTHONPATH=src:scripts PYTHONNOUSERSITE=1 \
  .venv/bin/python scripts/compare_p41z_jist_transfer_crop.py \
  --reference tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --atlas-gm build/p41z-case-normalization/iterations-0008-to-0011-surfaces/maseg-cortical-separation/Atlas_GM.nii.gz \
  --dura-dir build/p41z-case-normalization/iterations-0008-to-0011-surfaces/lesion-dura \
  --output-dir build/p41z-case-normalization/iteration-0011-java-crop-transfer
```

Exact inputs, per-map errors, crop bounds, Java class hashes, and selected
wrapper hash are in the [machine-readable comparison](p41z-jist-transfer-crop-comparison-2026-10-11.json).

## Next check

Evaluate source-compatible crop products through the same downstream fused
maps and surface-quality audit, then compare them with the recommended
physical-affine route. Keep the existing RAS-grid and ribbon acceptance
criteria unchanged.
