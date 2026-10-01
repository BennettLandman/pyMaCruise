# LesionTOADS transformed-template state handoff

The compressed NIfTI fixture is a synthetic 5 x 6 x 7 integer-label volume on
a positive-direction RAS grid (0.8 x 1.0 x 1.2 mm). The Java reference harness
passes the same byte volume to `updateTransformedTemplate` and records that the
preprocessing state retains the identical array reference and sets
`transformed=true`.

This fixture validates only the post-resampling state handoff. It does not
validate registration, label resampling, centroids, or segmentation.
