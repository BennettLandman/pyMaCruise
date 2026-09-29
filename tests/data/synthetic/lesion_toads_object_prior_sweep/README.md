# LesionTOADS OBJECT shape-prior sweep fixture

The deterministic `5×4×6` float32 phantom has three ordinary OBJECT classes,
positive data terms, six-face smoothing, asymmetric class smoothing priors,
off-diagonal similarities, and mixed DIRECT/PRIOR relations. Per-class
transformed shape values vary by coordinate; one class has zero support at an
interior coordinate to exercise the selected shape gate. The formula harness
matches the source in-place x/y/z scan, adds pairwise shape-prior energy,
converts and normalizes each voxel, then stores all interior memberships as
float32 bit patterns.

This fixture validates only the class-only OBJECT branch. It excludes lesion
and black-hole classes, auxiliary maps, registration, and full segmentation.
Regenerate with
`scripts/generate_java_lesion_toads_object_prior_sweep_reference.sh`.

Selected source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
