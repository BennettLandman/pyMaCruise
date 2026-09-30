# Selected LesionTOADS thinning-pass phantom

`thinning_pass_inputs.nii.gz` is an 11×11×11, 1-mm RAS synthetic volume
packed with 14 feature frames: segmentation, initial skeleton, memberships
for three classes, lesion and black-hole maps, transformed shape maps,
initial availability for each class, and an ordering field. Three x-axis
slabs make tissue classes with distinct volumes; class 1 is a lesion class.
Memberships favor the current class, so this fixture tests initialization,
source volume/minimum-volume accounting, max-heap evolution, per-class
cutoffs, and final state without changing labels through topology events.

The chosen `volumeRatio=0.4` makes class 0 use the ratio-derived minimum
above 100, while classes 1 and 2 use the source's 100-voxel minimum floor.
`java_reference.json` comes from the selected `initThinningLabels`,
`Numerics.max(100, volumeRatio*volume[k])`, and `propagateThinning` excerpts
plus the copied Java heap. Regenerate it with
`scripts/generate_java_lesion_toads_thinning_pass_reference.sh`; regenerate
the NIfTI with `scripts/generate_lesion_toads_thinning_pass_fixture.py`.

This validates one synthetic thinning pass with fixed memberships and shape
maps, not the upstream lesion-membership estimator or complete LesionTOADS.
