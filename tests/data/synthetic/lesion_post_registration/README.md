# LesionTOADS post-registration branch fixture

`branch_inputs_ras.nii.gz` stores the 16 fixed input frames from the P39zzzzs
statistics/membership phantom followed by synthetic segmentation, skeleton,
and ordering frames. Its affine is the source phantom affine and identifies
the voxel axes as RAS. The asymmetric segmentation contains a folded class
boundary; branch tests also derive a compact accepted/rejected topology queue
and an empty queue from these committed frames.

The companion Java reference pins the selected VandyCRUISE wrapper branch.
The component growth and classification maps are also checked against their
own Java-reference-backed tests, so this composition test focuses on call
order, float32 energy-ratio decisions, signed iteration-limit behavior, and
which downstream component runs.

To regenerate the phantom frames, run
`python tests/data/synthetic/lesion_post_registration/generate_fixture.py`.
The Java decision reference is emitted by
`tests/java_reference/generate_lesion_post_registration_reference.sh`.
