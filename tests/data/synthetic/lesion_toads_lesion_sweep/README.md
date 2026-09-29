# LesionTOADS LESION-class smoothing sweep fixture

This deterministic `5×4×6` float32 synthetic volume has two ordinary OBJECT
classes and one LESION class, a separate lesion-membership map and lesion data
term, coordinate-varying class data, mixed direct/prior relations, and one
interior zero-shape location. The Java formula harness follows the selected
VandyCRUISE x/y/z in-place scan, applies the selected class and LESION
six-face smoothing equations, converts memberships, and normalizes class and
auxiliary memberships at each voxel. Output float32 bit patterns are committed
for every interior voxel.

This is a bounded formula comparison only. Atlas prior additions, distance
modulation, black-hole/OUTMASK branches, registration, topology, and full
LesionTOADS segmentation remain excluded. Regenerate with
`scripts/generate_java_lesion_toads_lesion_sweep_reference.sh`.

Selected source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
