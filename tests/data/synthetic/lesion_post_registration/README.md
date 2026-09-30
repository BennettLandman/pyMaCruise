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

Pinned source provenance: `AlgorithmLesionToads.java` SHA-256
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`,
`LesionToadSegmentation.java` SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`, and
`LesionToadAtlas.java` SHA-256
`718ac87a670b1a9be1fcaf21d267be8d7150187c94003c62af0c8877534cb68e`. The
selected call-order excerpt and Java reference data are committed with this
fixture, so the regression test needs no external checkout.
