# TOADS relation-factor composition fixture

This fixture composes the selected `computeExactMemberships` support-tier
selection with its direct, neighbor, and prior factors on an asymmetric
synthetic four-slot label table. It records float32 factor bits for five
prior strengths and shape values that are negative, exactly zero, tiny
positive, and positive. Zero-shape prior classes receive zero weight as the
source's early membership exclusion requires; exact zero on direct or
neighbor support does not trigger that exclusion. There is no medical image
data.

Regenerate with
`scripts/generate_java_toads_relation_factor_composition_reference.sh`.
Selected-source SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884` for
`MultipleToadSegmentation.java`.
