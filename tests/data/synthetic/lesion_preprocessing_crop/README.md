# LesionTOADS preprocessing crop fixture

The fixture contains two normalized image channels and an asymmetric
three-label topology map, each as `.nii.gz` with an explicit RAS affine. The
image support and topology support differ so their union affects the inclusive
crop bounds. The active wrapper values are background label 1, background
ratio 0, and 10 zero-pad voxels. The axial RAS seed transform is identity.

`java_reference.json` stores the Java excerpt's crop bounds, padded dimensions,
cropped channel values as float32 bits, and transform update. Generate the
inputs with `generate_fixture.py` and the excerpt output with
`scripts/generate_lesion_preprocessing_crop_reference.sh`.

The selected source is VandyCRUISE `LesionToadPreprocess.java`, SHA-256
`72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08`.
This validates the fixed-normalized-input `findCroppingBoundaries()` and
`cropImages()` stages only. `normalizeImages()`, transform initialization,
topology resampling, uncropping, and registration remain outside this slice.
