# P39-to-plus transfer and crop fixture

This fixture validates the selected source-order bridge from the LesionTOADS
`dura removal inputs` output into the cropped images consumed by the
MaCRUISEplus pre-TGDM stages. It contains only synthetic, asymmetric RAS
NIfTI-1 images on a shared 17 x 17 x 17 grid; no subject data is used.

The test creates the selected P39ar stop state from the existing committed
folded LesionTOADS phantom, runs `lesion_toads_post_loop_outputs` with result
mode `dura removal inputs`, maps the returned GM, WM-fill, and WM-mask channels
onto the synthetic Seg reference, and crops them together with multi-atlas GM,
WM, mask, brain mask, and Seg. Multi-atlas GM determines bounds with the
source strict `> 0` support test and five-voxel border. All eight cropped
images have one affine and remain RAS. The CSF channel is intentionally not
returned: the selected `run_LesionToads_pipeline` call returns only GM,
WM-fill, and WM-mask; `run_LesionToads.m` writes CSF separately.

The Java reference harness covers the RAS identity-grid slice: source and
Seg grids have identical shape, affine, and RAS axis directions, so the
selected `ReorientVolume.solve(..., UNCHANGED, NEAREST_NEIGHBOR)` orientation
matrix is identity and voxel values are preserved before the target header is
assigned. The harness also calculates inclusive GM crop bounds, border
padding, and cropped values in source x/y/z order. It is a pinned Java
formula/reference harness, not execution of the dependency-bearing JIST
classes. The supplied source calls `ToadsCruiseSpaceToOriginalSpace` with
`ReorientVolume` and `CropBasedOnGM` with `CubicVolumeCropper`.

**Scope limitation:** the reusable Python transfer helper resamples in
physical coordinates onto the reference shape and affine. The selected JIST
method keeps source resolution/dimensions while orienting to target axis
directions and then assigns the reference header. This reference establishes
only identity-grid agreement; changed origins, shapes, or spacings are not
claimed equivalent. Full JIST/MIPAV runtime dependencies are not present in
the standalone test environment. This is not full plus or LesionTOADS
equivalence.

## Fixture files and SHA-256

| File | SHA-256 |
|---|---|
| `brain_mask.nii.gz` | `efcc9627df96d9ffab02eb43a221dead090048f01b70c9ed5f68828865ca0b3a` |
| `lesion_output_grid.nii.gz` | `901a055432781270df187d1bf3816e2003ddb26eff1532da05a92a4c8be0fed8` |
| `multi_atlas_gm.nii.gz` | `32a5ba2c080c45733157ef37508bd8188d25e29076e6d5058a2b1661f00212a1` |
| `multi_atlas_mask.nii.gz` | `e51fbafd8cb1c2519fc698d25b58c19efaab873a20244daf3c36f195e0278b57` |
| `multi_atlas_wm.nii.gz` | `9fc75ae8a7ee3e209d88634a16a155f06d631a6e0eaac9dd6dfcc258ba5828a1` |
| `segmentation_reference.nii.gz` | `83750fd1906b5a7b7d305dc991bde491a908de422bd3e220d20753a23973f0c7` |

Regenerate the NIfTI files with `python tests/data/synthetic/plus_lesion_transfer_crop/generate_fixture.py`.

## Source pins

| Selected source | SHA-256 |
|---|---|
| `ToadsCruiseSpaceToOriginalSpace.java` | `07e8b430944696568b0854c2352434f521e7939d99180be388cfaa2cf8a4656b` |
| `ReorientVolume.java` | `67fbf8b3f29e844cbec5e5994fecddf542da1e32d48222ee47312acb6886da24` |
| `CropBasedOnGM.java` | `4e41266b67a344cce53bc354b2302bf35486d1beb408e76a204d8c24d9fd1912` |
| `PlusLesionTransferCropRasIdentityReferenceHarness.java` | `a3fb7afa0365609e513d0cb480e169643637144e3010909bd0409befec266c94` |

The crop adapter reuses `pymacruise.cropping.find_crop_parameters` and
`crop_image`; the spatial handoff reuses
`pymacruise.spatial_transfer.transfer_to_reference_space`.
