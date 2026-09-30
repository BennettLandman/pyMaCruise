# P39x fixed-input repeated loop fixture

P39x reuses the folded asymmetric multi-channel RAS volume
[`iteration_inputs_ras.nii.gz`](../lesion_continuing_iteration/iteration_inputs_ras.nii.gz)
used by P39v/w. No external source checkout is read at test runtime. The
repeat-loop Java excerpt is in
[`LesionToadsFixedInputLoopReferenceHarness.java`](../../../java_reference/LesionToadsFixedInputLoopReferenceHarness.java)
and is generated into `java_reference.json` by
[`scripts/generate_lesion_fixed_input_loop_reference.sh`](../../../../scripts/generate_lesion_fixed_input_loop_reference.sh).

The Java excerpt pins the selected wrapper hash and repeated stop/growth
control over fixed energy histories. It does not duplicate the numerical
registration, statistics, membership, topology, or label-map kernels; those
retain their own selected-method Java comparisons. The image is synthetic,
RAS-oriented, and stored as `.nii.gz`.
