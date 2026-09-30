# TOADS membership relation-tier fixture

The fixture records selected Java exact-membership logic that assigns each
class at each voxel to direct support (closest-label slots 0/1), neighbor
support (slots 2/3), or prior support. It also records whether the source's
early zero assignment applies to a prior-tier class with an exactly zero
transformed shape. The small table includes missing slots, duplicate labels,
slot-priority conflicts, zero/tiny/negative/positive shape values, and all six
classes. There is no medical image data.

Regenerate with `scripts/generate_java_toads_relation_tiers_reference.sh`.
The selected source is `MultipleToadSegmentation.java`, SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
