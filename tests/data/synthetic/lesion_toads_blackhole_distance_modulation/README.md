# LesionTOADS black-hole distance-modulation reference

This synthetic selected-Java formula fixture contains fourteen scalar cases
for the `BLACKHOLE` branch's ventricle, gray-matter, and brainstem modulation.
It records lesion membership, black-hole auxiliary membership, and
`lesionRelation` output bits. In the selected source, a valid ventricle branch
scales lesion membership and relation but leaves black-hole membership
unchanged; a selected GM branch scales both auxiliary maps; a selected
brainstem distance then scales both maps independently.

The cases include ventricle-over-GM precedence, GM disabled/enabled, zero
distance with `spread`, `-1` sentinels, zero maxima, and independent brainstem
scaling. Regenerate with
`scripts/generate_java_lesion_toads_blackhole_distance_modulation_reference.sh`.
This fixture does not compute distance maps or normalize voxel memberships.
