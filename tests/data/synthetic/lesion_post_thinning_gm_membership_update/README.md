# Post-thinning statistics to GM-enabled membership update

This is a fixed-state continuation of the committed 17×17×17 folded RAS
phantom in `../lesion_initial_thinning_composition/initial_thinning_inputs_ras.nii.gz`.
It reuses the source-pinned `GMflag=false` first membership/thinning result and
the separate Java-pinned post-thinning statistics refresh. The next standard
13-class membership volume pass consumes those refreshed centroids, standard
deviations, and channel weights with `GMflag=true`.

Shape values, direct-relation flags, ventricle/GM/brainstem distances, and
initial atlas maps are fixed inputs from that fixture and the P39ak Java
reference. This slice does not generate those maps, re-run registration, or
claim a complete LesionTOADS iteration. The Java membership reference harness
consumes statistics output already compared bit-for-bit with its independent
P39al Java method harness, then applies the selected source-order volume
membership equations and writes all maps and objective bits.

This two-label segmentation leaves 11 OBJECT classes absent after thinning,
so their selected post-thinning standard deviations are zero. Java's
similarity arithmetic consequently produces 37,125 NaN class-energy voxels
and a NaN objective for this deliberately sparse phantom. The fixed
GM-enabled adapter preserves those computed IEEE results while still
rejecting caller-supplied nonfinite arrays; this is an explicit synthetic
edge-case outcome, not evidence that the profile yields a usable objective
for real inputs.

Regenerate the reference with
`python3 scripts/generate_lesion_post_thinning_membership_reference.py`.
