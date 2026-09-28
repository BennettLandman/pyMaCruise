# TOADS transformed shape-sampling reference

This fixture uses a synthetic 6×5×4 membership field on a template grid and
eight sample cases on a 7×6×5 input grid. It covers identity and fixed rigid
transform mapping, unequal image/template voxel sizes, source closest-edge
interpolation, slope recovery/clipping, and the class-specific defaults for
points outside their input-grid bounds. The shape values are generated from
`0.2 + 0.02*x + 0.025*y + 0.03*z + 0.005*x*y`.

The Java harness calls the copied selected VandyCRUISE source implementations
of `ParametricTransform.imageToTemplate`, `ParametricTransform.computeRotation`,
`ImageFunctions.linearClosestInterpolation`, and `Numerics.bounded`. The
source copies are hash-checked by the test suite. The harness applies the
source `MultipleToadDeformableAtlas.getRegularTransformedShape` class-bound
defaults and `shapeSlope` postprocessing around those methods; it does not
instantiate the full atlas class, run registration, or calculate optimized
class bounds. This is fixed-transform prior-sampling evidence only.

Reference source SHA-256 values:

- `ParametricTransform.java`: `22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac`
- `RotationMatrix.java`: `144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0`
- `ImageFunctions.java`: `be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892`
- `Numerics.java`: `13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb`

Regenerate with a JDK installed:

```sh
scripts/generate_java_toads_shape_sampling_reference.sh
```
