# LesionTOADS initial alignment bounds fixture

This fixture extends the existing `lesion_initial_atlas_alignment` synthetic
RAS subject through the selected post-pyramid
`LesionToadAtlas.computeTransformedShapeBoundingBox()` handoff. The already
validated initial-membership and rigid pyramid implementations are reused;
the new method computes only final transformed positive-support bounds.

Three committed RAS `.nii.gz` templates have class-dependent bounded support
and folded positive intensities. The Java reference receives the final rigid
transform from the Python composition, then uses the vendored selected
`ParametricTransform`, `RotationMatrix`, and `ImageFunctions` sources to scan
the transformed shape support. Raw transform bits, inclusive lower/upper
voxel coordinates, and support counts are compared. This directly validates
the source bounding-box scan and its interpolation boundary behavior.

The selected call sequence is
`registerShapesPyramid` → `computeTransformedShapeBoundingBox`. The initial
call order, seed, and membership bits remain covered by the adjacent fixture;
the numerical pyramid optimizer remains covered by its separate P39o
selected-Java comparisons. This does not claim whole-workflow or full
LesionTOADS equivalence.
