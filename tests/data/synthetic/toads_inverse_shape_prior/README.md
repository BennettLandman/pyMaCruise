# TOADS inverse shape-prior phantom

This fixture contains a deterministic 5×6×4×3 transformed-shape array stored
as a 1 mm RAS `.nii.gz` volume. Values use
`((3*x + 5*y + 7*z + 11*class) % 23) / 22`. The Java formula harness combines
it with shape priors `[0.2, 0.35, 0.45]`, the similarity matrix from the
committed P18ao source-formula fixture, and a prior coefficient of `0.15`.

`java_reference.json` records the source inverse-energy and half-weighted
objective-energy terms for every voxel/class, including zero border terms.
The selected `MultipleToadSegmentation.java` SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
Regenerate the image with `scripts/generate_toads_inverse_shape_prior_fixture.py`
and the formula output with
`scripts/generate_java_toads_inverse_shape_prior_reference.sh`.
