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

## Post-registration recrop sequence

`post_registration_recrop_java_reference.json` records a second crop lifecycle
from `LesionToadPreprocess`: first crop the committed two-channel image, undo
that crop and its preprocessing-transform translation, find bounds again from
the restored images and topology, and crop a second time. The source-excerpt
harness is `tests/java_reference/LesionPostRegistrationRecropReferenceHarness.java`.
It compares inclusive bounds, crop dimensions, float32 transform bits, and the
SHA-256 digest of every float32 crop voxel. Array bytes use channel/x/y/z order
and little-endian float32. Regenerate it with
`tests/java_reference/generate_lesion_post_registration_recrop_reference.sh`.

The copied crop, uncrop, and boundary method behavior is pinned to
`LesionToadPreprocess.java` SHA-256
`72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08`
(selected source lines 344–479); coordinate mapping is lines 653–660.
Supporting rotation and floor/ceil excerpts use `RotationMatrix.java`
`144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0` and
`Numerics.java`
`13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb`.
This is a selected-method harness, not execution inside the full MIPAV/JIST
application. The Python test also asserts that the optimized atlas transform
is independent of the preprocessing transform used for recrop bounds.
