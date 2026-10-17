# AtlasEM Dura-input NIfTI handoff — 2026-10-16

`run_dura_from_atlas_em_products` accepts the native-grid compressed NIfTI
created by the AtlasEM execution adapter. It validates a three-channel RAS
product (`GM membership`, `WM membership`, `filled WM mask`), checks its grid
against the original T1, and passes the channels in the selected JIST input
order to the existing Dura implementation with the wrapper's V3 parameters.
`run_atlas_em_prelesion_dura` now uses the same function after generating its
products from memberships.

The full synthetic descriptor case in
[`atlasem-to-dura-full-synthetic-2026-10-16.md`](atlasem-to-dura-full-synthetic-2026-10-16.md)
passes the restored `atlasem_dura_inputs.nii.gz` image to this interface and
completes one Dura outer pass. The previous independently generated tissue
phantom also gives the same Dura brain-mask output through the product-image
interface as through the membership-map convenience API. Focused Dura/CLI
tests pass **2 tests, 3 expected RAS warnings, in 0.90 s**; the broader
AtlasEM/registration/Dura/morphology group passes **44 tests, 3 expected RAS
warnings, in 6.35 s** on test commit
`c5abb414fe504a0f280ec2806ad8f32c204b66cb`.

This validates the synthetic file-level handoff and parameter wiring, not
selected-atlas or human pre-LesionTOADS parity.
