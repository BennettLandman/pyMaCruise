# LesionTOADS OBJECT in-place sweep fixture

The deterministic `5×4×6` float32 fixture has three ordinary OBJECT classes,
direct relations, positive shape support, positive fixed data terms, positive
class smoothing priors, and six-face smoothing. A Java formula harness follows
the selected `computeMemberships` x/y/z loop, writes each voxel immediately,
normalizes its class memberships before moving on, and records every interior
voxel's final class memberships as float32 bit patterns. The coordinate-distinct
initial memberships ensure that scan order affects later neighborhoods.

This isolates the source in-place update and normalization order only. It does
not include atlas priors, non-OBJECT relations/types, lesion auxiliaries,
topology competition, registration, or complete LesionTOADS segmentation.
Regenerate with `scripts/generate_java_lesion_toads_object_sweep_reference.sh`.

Selected source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
