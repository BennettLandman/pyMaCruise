# LesionTOADS lesion distance-modulation reference

The selected Java formula fixture has twelve scalar cases for modulating
lesion membership and its relation factor by already-computed ventricle, gray
matter, and brainstem distances. It exercises ventricle-over-GM precedence,
the GM enable flag, zero-distance `spread`, `-1` sentinels, zero maxima, and
the independent brainstem multiplier.

Regenerate with
`scripts/generate_java_lesion_toads_distance_modulation_reference.sh`. This
does not compute distance maps, modulate black-hole membership, or normalize
voxelwise memberships.
