# TOADS relation-map diffusion fixture

Five synthetic Java references cover a uniform segmentation, planar interface,
and folded six-class label field. Cases include zero support, zero limit, and
positive membership support with iteration-cap and fixed-point behavior. Each
reference records float32 membership inputs, the initial relation count, the
number of Java loop passes, and the final map. No medical images are used.

The harness composes selected excerpts of `buildRelationList` and
`computeRelations` from
`VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`,
SHA-256 `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
Regenerate with `scripts/generate_java_toads_relation_diffusion_reference.sh`.
The tested path is exact-membership mode; it does not cover approximate mode,
topology evolution, or full TOADS segmentation.
