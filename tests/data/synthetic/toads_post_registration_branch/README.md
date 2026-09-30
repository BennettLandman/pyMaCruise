# TOADS post-registration branch fixture

This fixed-state fixture exercises only `AlgorithmMultipleToads.runAlgorithm`
after the selected default-atlas membership result and rigid registration are
available. It is synthetic; it is not an imaging-data validation.

`membership_maps.nii.gz` and its objective are copied from case 5 of the
selected Java-reference default 11-class membership sweep. The fixture tests
verify the membership bits against that reference. `segmentation.nii.gz` is an
asymmetric folded, three-label XYZ volume. `skeleton.nii.gz` has one
deterministic interior growth site. The ordering, class volumes, boundary
distances, and prior output labels are fixed branch inputs.

All NIfTI files use a positive-determinant RAS affine with anisotropic voxel
spacing. The generated registration transform is a fixed synthetic handoff
value and the branch must carry it through without mutation. The separate
integration test also feeds the real `ToadsDefaultRegistrationHandoff`
produced by P18czh through this branch adapter.

`java_reference.json` is emitted by
`scripts/generate_java_toads_post_registration_branch_reference.sh`. The
harness pins the selected `AlgorithmMultipleToads.java` SHA-256 and records the
float32 energy-ratio bits, stop/grow flags, and source-ordered branch stage
dispatch for convergence, iteration ceiling, and both `doGrowing` settings.
The growth, classification, and relation-state array kernels are reused from
their separately pinned Java-reference tests; this harness validates their
outer-loop branch selection and call order.
