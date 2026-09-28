# Folded TOADS shape support phantom

This synthetic atlas shape prior is a 39×35×31 float32 RAS-oriented `.nii.gz`
volume with 1 mm isotropic voxel spacing. It contains a smooth, folded,
cortex-like ribbon around an ellipsoidal shell. Its deterministic generator is
`scripts/generate_folded_toad_shape_prior.py`; it uses only analytic functions
and does not represent a human image.

`java_reference.json` records the positive-support count and exact inclusive
lower/largest-positive-index bounds found by the selected Java
`computeTransformedShapeBoundingBox` scan excerpt for identity, translated,
and rotated-plus-translated rigid transforms on a 43×39×35 input grid. The
harness calls the selected Java `ParametricTransform` and
`ImageFunctions.linearClosestInterpolation` APIs. The corresponding Python
check loads the committed NIfTI prior and recomputes each bound. This validates
bound refresh for these synthetic transforms; it does not validate optimizer
alignment or TOADS segmentation.

Regenerate the Java formula reference with
`scripts/generate_java_toads_folded_shape_bounds_reference.sh`.
