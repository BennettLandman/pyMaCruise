# TOADS boundary-distance fixture

Selected Java references exercise `computeBoundaryDistances` on a planar and
a folded synthetic categorical volume. They include fresh and prepopulated
signed-byte distance arrays and limits 0, 1, and 3. The fixture intentionally
records that the selected source seeds both class sides at a six-face
interface but does not spread the values: its inner `min(0, positive + 1)`
remains zero, and its loop `change` flag is cleared without being restored.
There is no medical image data.

Regenerate with `scripts/generate_java_toads_boundary_distances_reference.sh`.
The source is
`VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`,
SHA-256 `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
The comparison covers this method only, not TOADS distance-weighted
membership behavior or complete segmentation.
