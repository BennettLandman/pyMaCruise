# LesionTOADS transformNewImage fixture

This fixture validates only the selected, dormant `LesionToadPreprocess.transformNewImage()` fixed-transform sampler. It maps a template-grid intensity image to an inclusive target-image crop, initializes the padded frame to the caller's background value, computes the Java inverse-Rodrigues template coordinates, and uses `ImageFunctions.linearInterpolation` with the selected background behavior at invalid interpolation coordinates. The selected `AlgorithmLesionToads` wrapper contains no call site for this utility; this is not evidence of active workflow behavior or full LesionTOADS equivalence.

`template_image.nii.gz` is a synthetic `(9, 10, 11)` float32 volume with RAS axes, voxel sizes `(0.9, 1.0, 1.1)` mm, and a nonzero RAS origin. Its asymmetric values use a deterministic integer polynomial divided by 16 so all source samples are exactly representable as float32. The harness applies image shape `(7, 8, 9)`, inclusive bounds `(1,5,1,6,1,7)`, one-voxel padding, background `-2.5`, image spacing `(1.1,1.3,0.9)`, template spacing `(0.9,1.0,1.1)`, and transform `(0.02,-0.03,0.01,0.3,-0.15,0.4)`. It records two complete float32 bit-pattern outputs: identity direction signs `(1,1,1)` and a second case with negative first-axis direction `(-1,1,1)`. The NIfTI fixture remains RAS in both comparisons; the signed factor tests the selected source sampling geometry.

Selected-source SHA-256 values:

- `LesionToadPreprocess.java`: `72a1a47dae9c448022707ee68917eb5ca70d64b763f38582fa31304463732c08`
- `ImageFunctions.java`: `be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892`
- `RotationMatrix.java`: `144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0`

Fixture hashes: harness `6676ff3cc7393e040e4affa9bd3e9e009701395e092fdfeccf605f747e6ebca0`, Java reference `ab8eaedcf39dbec3e686509827eeb5773a91a461320f5c903a39f144d49e6b00`, and NIfTI input `23d50b18fc1d0ae4ecd92683791f7087af15d208f2c3c676bdd68fe0b5054caa` (SHA-256).

The reference now also includes a second translated output that drives valid
crop voxels outside the template's interpolation domain. This distinguishes
the interpolation background path from the crop's initially padded faces.

Regenerate the in-repository reference and NIfTI fixture with:

```sh
mkdir -p tmp/build/lesion-transform-new-image
javac -d tmp/build/lesion-transform-new-image \
  tests/java_reference/LesionToadsTransformNewImageReferenceHarness.java
java -cp tmp/build/lesion-transform-new-image \
  LesionToadsTransformNewImageReferenceHarness \
  > tests/data/synthetic/lesion_transform_new_image/java_reference.json
PYTHONPATH=src .venv/bin/python \
  tests/data/synthetic/lesion_transform_new_image/generate_fixtures.py
```

The fixture generator writes only `template_image.nii.gz` in this directory.
