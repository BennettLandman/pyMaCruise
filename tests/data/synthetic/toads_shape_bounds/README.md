# TOADS transformed shape support bounds

`java_reference.json` stores five synthetic support scans on a 6×5×4 input
grid using a 6×5×4 prior. It covers identity and rigid transforms at unequal
voxel resolutions, translation, support touching all image faces, and an
atlas prior with no positive samples. The positive box prior has negative
values outside its support so trilinear interpolation exercises the `> 0`
support threshold.

The harness calls the copied selected Java
`ParametricTransform.precomputeImageToTemplateMatrix` and
`ImageFunctions.linearClosestInterpolation` methods, then applies the
`MultipleToadDeformableAtlas.computeTransformedShapeBoundingBox` scan excerpt.
The atlas source SHA-256 is
`e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`; the
transform source SHA-256 is
`22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac`.

The scan preserves two source details: it initializes empty support to
`lower=image_dimensions, upper=(0,0,0)`, and it stores the greatest support
index without adding one. Later source sampling treats that maximum as an
exclusive bound, excluding the maximum-support plane. The fixture verifies
the bounds and this lookup behavior. It does not port atlas registration,
class-bound optimization across iterative transforms, or physical-world
origin reconciliation. All data are synthetic numeric arrays, not MRI.

Regenerate with a JDK installed:

```sh
scripts/generate_java_toads_shape_bounds_reference.sh
```
