# LesionTOADS constrained-centroid phantom

`centroid_inputs.nii.gz` is a deterministic 5×4×3, 1-mm RAS synthetic volume
with nine float32 frames: hard class labels, four membership maps, two image
channels, and two correction fields. Labels 0–2 occur; class 3 is empty.
The patterned labels/intensities/memberships exercise source x/y/z accumulation
order, and the data means exceed the initial centroid bands for classes 0 and
1 so the source upper clamps are used. The empty class tests the source's zero
mean followed by its lower clamp.

Selected-method Java harnesses provide raw float32 references for both
`computeConstrainedCentroids()` and `computeWeightedCentroids(factor)`. The
former tests field/no-field behavior and 45% bounds. The latter tests factors
0, 0.35, and 1 with and without `useField`, including the source behavior
that an empty class becomes zero even at factor 1. Python matches all bits,
including field multiplication before accumulation and the x/y/z traversal.
Regenerate the NIfTI with `scripts/generate_lesion_centroid_fixture.py`; the
Java references are regenerated with
`scripts/generate_java_lesion_constrained_centroids_reference.sh` and
`scripts/generate_java_lesion_weighted_centroids_reference.sh`.

These validate two centroid update methods on fixed labels and memberships.
They do not estimate lesion/black-hole centroids, variances, channel weights,
membership energies, or a full LesionTOADS iteration.
