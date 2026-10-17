# AtlasEM native-grid RAS products — 2026-10-16

`restore_atlas_em_native_nifti_products` packs the final crop's selected
AtlasEM memberships, three Dura-input channels, internal mask, and tissue-index
segmentation; restores them to the native T1 dimensions using the source-shaped
uncrop; and assigns the normalized reference T1's affine/forms. The companion
`save_atlas_em_native_nifti_products` writes four `.nii.gz` files:

- `atlasem_memberships.nii.gz` — class-last float32 membership maps.
- `atlasem_segmentation.nii.gz` — signed int8 tissue indices, including the
  source's `-1` empty sentinel.
- `atlasem_internal_mask.nii.gz` — uint8 AtlasEM statistical mask.
- `atlasem_dura_inputs.nii.gz` — float32 GM, WM, and filled WM-mask channels.

The synthetic descriptor integration checks that each file is written, uses
`.nii.gz`, returns to the full native grid, and carries an RAS affine equal to
the supplied RAS T1. The broader focused AtlasEM/registration/Dura/morphology
group passed **43 tests, 3 expected RAS orientation warnings, in 6.38 s** on
implementation commit `4b123555c266be02ec8910b47774a7e6cc91bfd4`.

This covers output geometry, dtype, filenames, and compressed NIfTI persistence
for the synthetic execution. It does not validate data against the missing
selected `ibsr8obj` atlas or run the Dura evolution on those generated maps.
