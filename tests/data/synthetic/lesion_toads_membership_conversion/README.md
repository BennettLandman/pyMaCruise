# LesionTOADS class membership conversion reference

The selected Java fixture contains eleven scalar formula cases for converting
a completed class data energy to membership. It exercises values below, at,
and just above `ZERO`, the finite `INF` sentinel, direct and prior relations,
and `shapeFactor = Numerics.min(1, 2*shape)` below, at, and above saturation.

Regenerate with
`scripts/generate_java_lesion_toads_membership_conversion_reference.sh`.
This validates one class conversion only, before lesion auxiliary terms,
distance/relation modulation, denominator accumulation, or normalization.
