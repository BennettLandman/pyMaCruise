# LesionTOADS thinning outer-loop phantom

`queue_run_inputs.nii.gz` is a deterministic 5×5×5, 1-mm RAS synthetic
volume. Fourteen float32 feature frames store segmentation, skeleton, three
membership maps, lesion and black-hole maps, three transformed-shape maps,
three availability maps, and ordering values. The test constructs the
selected max heap and scalar class-volume thresholds explicitly.

The event run propagates a front through a small all-class-0 interior,
discards one stale skeleton item and one class at its minimum volume, and
includes an initial item exactly at `-limit` plus a lower-priority item. The
Java excerpt demonstrates the selected outer-loop behavior: because the loop
compares the previous priority, it pops and processes one crossing item
before it stops, leaving the lower-priority heap entry untouched. Build the expected JSON with
`scripts/generate_java_lesion_toads_thinning_queue_run_reference.sh`; rebuild
the NIfTI with `scripts/generate_lesion_toads_thinning_queue_run_fixture.py`.

This validates repeated queue evolution from supplied initial state. It does
not combine frontier initialization or the source `volumeRatio` minimum
calculation into a complete workflow pass.
