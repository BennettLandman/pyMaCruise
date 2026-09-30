# Composed LesionTOADS thinning-event phantom

`event_inputs.nii.gz` is a deterministic 5×5×5, 1-mm RAS synthetic volume
packed along the fourth NIfTI axis as segmentation, skeleton, three class
membership maps, lesion membership, black-hole membership, three transformed
shape maps, three per-class availability maps, and ordering values. Its 3×3×3
interior neighborhood is the topology-accepted direct-change case from the
P39zzi reference; all six face neighbors are available for class 0.

The fixture compares a max-heap with a selected center item and an unrelated
remaining seed. It checks the composed selected Java sequence: pop, target
selection, topology test, center state and volume update, lesion-aware
neighbor scheduling, and final queue pop order. Rebuild with
`scripts/generate_lesion_toads_composed_thinning_event_fixture.py` and
`scripts/generate_java_lesion_toads_composed_thinning_event_reference.sh`.
This is a single queue event; it does not exercise frontier initialization,
repeated thinning, or the outer limit loop.
