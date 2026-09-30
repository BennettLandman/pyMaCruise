# LesionTOADS thinning-neighbor scheduling phantom

`schedule_inputs.nii.gz` is a
deterministic 5×5×5, 1-mm RAS fixture packed as seven sets of six float32
feature volumes along the fourth NIfTI axis: memberships for classes 0–2,
lesion membership, black-hole membership, and the availability map for the
processing class. The test reconstructs the class axis and keeps all data
synthetic.

Cases cover six available faces for an ordinary class, lesion `useLesionMax`
and lesion-sum scoring, equal-ratio priorities exactly at and just above the strict `-limit`
cutoff, class volume equal to or below its minimum, and a pre-unavailable face.
The Java reference runs the selected neighbor loop and copied
`BinaryHeap4D`, recording every insertion, final availability, and heap pop
order (including equal-priority ties). Regenerate the image and JSON with
`scripts/generate_lesion_toads_thinning_schedule_fixture.py` and
`scripts/generate_java_lesion_toads_thinning_schedule_reference.sh`.

This is a scheduler slice, not a complete LesionTOADS thinning pass.
