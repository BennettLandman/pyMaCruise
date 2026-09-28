# TOADS shape-weighted membership comparison

This composition reuses committed RAS `.nii.gz` inputs from
`toads_direct_membership/` (data energy and initial memberships) and
`toads_shape_relations/` (shape priors and four-slot closest labels), avoiding
duplicate image fixtures. The tests load all four inputs and check their RAS
1 mm geometry.

`java_reference.json` contains the direct/neighbor/prior relation weighting,
prior-only zero-shape suppression, `ZERO=1e-10f` reciprocal rule, and
class-order normalization from a Java formula harness. The source hashes and
generator formulas are included in the JSON and the neighboring fixture
READMEs. Regenerate the reference with
`scripts/generate_java_toads_shape_weighted_membership_reference.sh`.

This is a per-voxel composition with the source's smoothing coefficient set
to zero; it does not execute the asynchronous volume sweep.
