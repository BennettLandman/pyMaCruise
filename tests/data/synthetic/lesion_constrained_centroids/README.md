# LesionTOADS constrained-centroid phantom

`centroid_inputs.nii.gz` is a deterministic 5×4×3, 1-mm RAS synthetic volume
with nine float32 frames: hard class labels, four membership maps, two image
channels, and two correction fields. Labels 0–2 occur; class 3 is empty.
The patterned labels/intensities/memberships exercise source x/y/z accumulation
order, and the data means exceed the initial centroid bands for classes 0 and
1 so the source upper clamps are used. The empty class tests the source's zero
mean followed by its lower clamp.

The selected-method Java harness calculates the weighted means and 45% bounds
with and without `useField`; it emits raw float32 centroid bits. Python matches
all centroid bits, including field multiplication before residual
accumulation, all exercised clamps, and empty-class behavior. The fixture and
reference can be regenerated with
`scripts/generate_lesion_centroid_fixture.py` and
`scripts/generate_java_lesion_constrained_centroids_reference.sh`.

This validates the constrained centroid updater on fixed labels and
memberships. It does not estimate lesion/black-hole centroids, variances,
channel weights, membership energies, or a full LesionTOADS iteration.
