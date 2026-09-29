# LesionTOADS rigid shape-prior sampling reference

The 7×6×5 synthetic input grid samples two 6×5×4 atlas shape maps. The
committed Java reference covers identity and nonzero rigid mappings, an
anisotropic image grid, fractional mapped coordinates, input bounds, and
coordinates outside the shape-map grid. It checks both the first structure
(class index 0, outside value 1) and other structures (outside value 0).

This follows the selected VandyCRUISE
`LesionToadAtlas.getRegularTransformedShape` wrapper behavior around its
`ParametricTransform.imageToTemplate` and `ImageFunctions.linearInterpolation`
calls. Unlike the TOADS sampler, it does not use closest-edge interpolation or
`shapeSlope` recovery. The source method returns its class-specific default
whenever any mapped coordinate is outside `[0, dimension-2]`.

Selected source `LesionToadAtlas.java` SHA-256:
`718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e`.
The Java formula harness is hash-checked by `tests/test_lesion_shapes.py`; the
shared Java transform/interpolation source-copy hashes are checked by the
existing TOADS shape-sampling tests. The harness mirrors this selected method's
small formula wrapper and does not instantiate the full atlas class, perform
registration, or exercise multiple/deformable transforms.

Regenerate from the repository root with
`scripts/generate_java_lesion_shape_sampling_reference.sh` (requires a JDK).
No clinical image data is used.
