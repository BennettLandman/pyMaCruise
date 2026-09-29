# LesionTOADS registration pyramid sampling fixture

The Java harness compares the selected rigid `ParametricTransform` scaled
image-to-template map with `ImageFunctions.linearInterpolation` and its three
spatial derivative methods. It uses an odd 21×19×17 synthetic input and a
24³ synthetic template, giving floor-truncated coarse grids for scales 2, 4,
and 8. The shape prior is a deterministic affine/bilinear ramp. A fractional
interior transform and large negative/positive translation cases cover
interpolation and constant-zero exterior behavior at each scale.

The full-resolution image center is `(10.5, 9.5, 8.5)` and the synthetic
template center is `(8, 8, 8)`; both remain in original-resolution units while
the arrays and sample coordinates are coarse. Selected source Java files are
compiled from checked-in copies by
`scripts/generate_java_lesion_scaled_shape_sampling_reference.sh`.
`tests/test_lesion_shapes.py` hash-checks both the harness and fixture.

This is formula-level registration pyramid evidence. It does not instantiate
`LesionToadAtlas`, calculate optimized class bounds, or run registration
optimization or lesion segmentation. No image assets or clinical data are
used.
