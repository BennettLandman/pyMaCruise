# TOADS OUTLIER voxel-normalization reference

Three selected-Java scalar cases validate TOADS's final denominator,
class-membership normalization, repeated OUTLIER auxiliary-map division, and
per-class objective contribution. The reference uses formula-defined
memberships for OBJECT, MASK, and OUTLIER classes; one case has two OUTLIER
classes to preserve source repeated-update semantics, and one has zero
denominator. The selected source hash is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
No image data are used. Regenerate with
`scripts/generate_java_toads_outlier_normalization_reference.sh`.
