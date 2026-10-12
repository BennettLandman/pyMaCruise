# P41z JIST native-transfer geometry audit — 2026-10-11

## Result

I compared the iteration-11 synthetic Dura maps against the selected JIST
`ReorientVolume.solve(UNCHANGED, NEAREST_NEIGHBOR)` using the saved T1 as the
reference. The Dura maps are RAS volumes of **71×89×67**; the reference T1 is
RAS **81×101×77**.

The Java output retains the source map extents (**71×89×67**) and its values
are exactly the source voxels for all four channels after converting MIPAV's
A2P storage order back to NIfTI RAS. The selected wrapper then assigns the
reference header. The Python `transfer_to_reference_space` adapter instead
uses the source and reference physical affines to sample a full **81×101×77**
volume. Its output shape therefore differs from the selected Java behavior on
this cropped case. This is why the audit records geometry and array behavior
separately; it does not claim voxelwise equivalence between these outputs.

This is a selected-wrapper geometry difference, not proof that either result
explains P41z's failed ribbon. The physical-affine Python result preserves the
crop's physical position, while the selected JIST call preserves the source
extents and subsequently assigns the reference header. The source's
`CropBasedOnGM` consumer must be evaluated before choosing a compatibility
policy. Keep physical-RAS transfer as the recommended behavior pending that
check; a future legacy mode could reproduce the selected wrapper when users
need matching prior outputs. The already validated same-grid human comparison
remains bitwise exact and is unaffected.

## Reproduction

The reusable audit compiles/runs the selected Java harness and reports Java
output extents, source-array preservation, and Python physical-transfer
grid. It requires the development MIPAV/JIST runtime artifacts and does not
add them as Python runtime dependencies:

```sh
PYTHONPATH=src:scripts PYTHONNOUSERSITE=1   .venv/bin/python scripts/audit_jist_reorient_native_transfer.py   --reference tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz   --map csf=build/p41z-case-normalization/iterations-0008-to-0011-surfaces/lesion-dura/lesion_toads_csf.nii.gz   --map gm=build/p41z-case-normalization/iterations-0008-to-0011-surfaces/lesion-dura/lesion_toads_gm.nii.gz   --map wm_fill=build/p41z-case-normalization/iterations-0008-to-0011-surfaces/lesion-dura/lesion_toads_wm_fill.nii.gz   --map wm_mask=build/p41z-case-normalization/iterations-0008-to-0011-surfaces/lesion-dura/lesion_toads_wm_mask.nii.gz   --output-dir build/p41z-case-normalization/iteration-0011-native-transfer-geometry-audit
```

The recorded input hashes, Java source/harness hashes, and per-channel
geometry are in the [machine-readable audit](p41z-jist-native-transfer-geometry-2026-10-11.json).
This is synthetic selected-wrapper evidence only, not whole-workflow source
equivalence or real-MRI validation.

## Next check

Run the selected Java crop consumer with its actual `Atlas_GM` reference and
compare its cropped maps against the Python full-reference transfer followed
by MaSeg-GM crop. That determines whether the source wrapper's retained extent
is consumed as historical behavior or exposes another workflow discrepancy.
