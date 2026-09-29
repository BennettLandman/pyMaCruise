# LesionTOADS distance modulation volume fixture

The committed selected-Java formula reference extends the OBJECT+LESION
shape-prior sweep to synthetic ventricle, GM, and brainstem distance maps.
The coordinate-dependent signed-byte maps exercise distance-zero spread,
ventricle precedence over GM, GM fallback with the flag enabled/disabled, and
independent brainstem modulation. Direct/prior relations and a zero-shape
LESION voxel remain in the phantom. Every interior class and lesion membership
is stored as a float32 bit pattern.

Regenerate with `scripts/generate_java_lesion_toads_lesion_distance_sweep_reference.sh`.
This fixture does not cover BLACKHOLE-specific auxiliary distance behavior,
other lesion classes, registration, or complete segmentation.

Selected source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
