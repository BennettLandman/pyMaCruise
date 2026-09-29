# LesionTOADS six-face neighborhood fixture

This committed reference uses deterministic coordinate-distinct float32
class, lesion, and black-hole membership volumes with shape `5×4×6` and three
class channels. The selected VandyCRUISE
`LesionToadSegmentation.computeMemberships` source loops over interior
coordinates in x/y/z order and reads neighbors in `X+, X-, Y+, Y-, Z+, Z-`
order. The Java formula harness records each selected voxel and its face
values as float32 bit patterns. It does not instantiate the segmentation
engine, compare an update, or establish segmentation equivalence.

Regenerate with:

```sh
scripts/generate_java_lesion_toads_six_face_neighborhood_reference.sh
```

Source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
