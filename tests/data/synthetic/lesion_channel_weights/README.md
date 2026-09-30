# LesionTOADS channel-weight phantom

`channel_weight_inputs.nii.gz` is a deterministic 5×4×3, 1-mm RAS volume
with six membership maps, two image channels, two field maps, and lesion,
black-hole, and outlier membership maps. The selected-method Java reference
uses OBJECT, OUTLIER, LESION, and BLACKHOLE contributors plus MASK and OUTMASK
classes, which this source method ignores.

The Java excerpt and Python implementation match raw float32 channel-weight
bits with and without field correction. The fixture tests the selected
outlier-distance term and the lesion-related linear residual formula (which
is not squared in the source). It also tests the non-diagonal uniform-weight
path and an empty-contribution case: the selected `INF=1e10` values normalize
to equal channel weights. Inputs are not mutated.

Regenerate the input with `scripts/generate_lesion_channel_weight_fixture.py`
and the reference with
`scripts/generate_java_lesion_channel_weights_reference.sh`. This validates
the channel-weight update only, not the upstream statistics, membership
iteration, or complete LesionTOADS workflow.
