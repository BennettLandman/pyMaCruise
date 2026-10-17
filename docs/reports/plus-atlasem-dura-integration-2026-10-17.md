# MaCRUISEplus AtlasEM/Dura handoff integration — 2026-10-17

## Change

The raw-input restartable PLUS runner accepts an optional distinct
`--atlasem-descriptor`. When supplied, it inserts a restartable
AtlasEM/Dura stage after ANTsPy FLAIR-to-T1 registration, N4 correction, and
the user-supplied brain mask, and before LesionTOADS intensity normalization
and crop. The AtlasEM algorithm runs on prepared T1; its selected Dura output
strips T1 and supplies the mask applied to registered FLAIR. The resulting
T1/FLAIR channel-last `.nii.gz` on the native RAS T1 grid feeds the existing
LesionTOADS crop.

The AtlasEM descriptor is separate from the LesionTOADS atlas descriptor.
Every referenced AtlasEM asset is validated and hashed in the stage
manifest; child outputs are hash-verified on resume. Omitting
`--atlasem-descriptor` retains the earlier LesionTOADS-only path. The Python
implementation contains the stage; execution on the selected
`ibsr8obj` descriptor is still blocked because that resource set has not
been staged. No generated test atlas is treated as a production substitute.

## Synthetic validation

A generated three-class AtlasEM descriptor runs one iteration at each
registration/EM/Dura limit. The test checks native-grid RAS `.nii.gz`
outputs, Dura-stripped T1 and Dura-masked FLAIR channel ordering, zeroed
outside-mask/negative FLAIR, stage manifest order, and same-input resume.
A parent-runner test confirms that, when enabled, the AtlasEM/Dura stage
completes before lesion crop and that its conditioned modalities are the
inputs passed to crop. The CLI parser accepts the optional descriptor.

Focused handoff, placement, parser, and legacy-route regression: **10 passed
in 66.29 s**, with four expected orientation warnings.

Full pytest on implementation commit
`51e87d3b15b6792313edb1a3b1783880b276778f`: **1,673 passed, 4 skipped,
0 failed, 21 expected orientation warnings in 665.27 s**. Peak sampled RSS
was **15.388 GB** under the 16 GB watchdog; it did not terminate the run. See
the [full-suite report](full-pytest-2026-10-17-plus-atlasem-dura-integration.md).

This validates synthetic stage composition and compatibility with the
existing runner. It does not validate the selected descriptor/payload,
human AtlasEM convergence, or complete PLUS source equivalence.
