# LesionTOADS smoothing/prior membership integration

Two selected-Java synthetic voxel cases compose the six-face smoothing terms
and shape/similarity prior terms before direct/prior inverse membership,
lesion or black-hole distance modulation, denominator normalization, and
local energy. One case exercises LESION with a direct relation and ventricle
plus brainstem distances; the other exercises BLACKHOLE with a prior relation,
the extra black-hole lesion-class prior, a zero GM distance, and brainstem.

Shape/similarity arrays and neighbor stencils are precomputed fixture inputs.
This is not the full volume update or construction of those spatial fields.
Regenerate with
`scripts/generate_java_lesion_toads_smoothing_prior_membership_reference.sh`.
