# TOADS classification transfer fixture

This fixture transcribes the selected `MultipleToadSegmentation.updateClassificationFromLabels`
loops on an asymmetric synthetic 3×2×4 categorical volume. The reference
stores the template label bytes after lookup and Java's changed-voxel count,
covering every template index, initial zero output, unchanged and changed
bytes, and boundary voxels. It contains no medical images.

Regenerate using `scripts/generate_java_toads_classification_labels_reference.sh`.
The selected source is
`VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`,
SHA-256 `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
The reference is a method-formula comparison, not whole TOADS segmentation
equivalence.
