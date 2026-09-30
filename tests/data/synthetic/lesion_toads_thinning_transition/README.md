# LesionTOADS thinning-event phantom

`transition_inputs.nii.gz` is a deterministic 3×3×3, 1-mm RAS synthetic
neighborhood pack. It stores seven cases, each as 11 float32 feature volumes
along the fourth NIfTI axis: segmentation, skeleton, ordinary membership for
classes 0–2, lesion membership, black-hole membership, transformed shape for
classes 0–2, and the popped queue priority. The test treats integer-valued
segmentation and skeleton feature frames as Java short labels. Initial class
volumes, minimums, availability, object types, processing classes, and
`useLesionMax` are deterministic test controls documented in the test.

The cases exercise direct topology accept/reject, the selected shape-prior
fallback, a target with no connected face neighbor, a voxel already assigned
its best class, the strict `volume[k] > minVolume[k]` gate, and lesion-class
volume reduction with `useLesionMax=true`.

`java_reference.json` is generated from the selected
`LesionToadSegmentation.propagateThinning` event excerpt and the already
validated LUT-backed topology excerpt. Rebuild it with
`scripts/generate_java_lesion_toads_thinning_transition_reference.sh`. The
reference is a local event comparison, not a complete thinning pass or a
claim of LesionTOADS equivalence.
