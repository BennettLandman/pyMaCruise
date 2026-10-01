# LesionTOADS rigid crop fixture

This companion case to `lesion_preprocessing_crop` holds the same folded,
asymmetric normalized two-channel image and three-label topology in a true RAS
NIfTI grid, but supplies a nonzero Rodrigues rotation and translation. It
therefore checks that topology support is mapped into image coordinates before
floor/ceil bounds are chosen and that crop translation is rotated when the
inclusive image crop is applied.

The Java reference is a selected-method excerpt for
`LesionToadPreprocess.findCroppingBoundaries()` and `cropImages()` at source
SHA-256 `72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08`.
The Java excerpt SHA-256 is
`43c1846cb48d2dec52dad9de1c85c805b7666bfeba5213335b950aaaaebda4e0`.
Generate inputs with `generate_fixture.py`, then generate the reference with
`tests/java_reference/generate_lesion_preprocessing_crop_rigid_reference.sh`.
This validates only the fixed-input crop methods, not normalization,
registration, topology resampling, or segmentation.
