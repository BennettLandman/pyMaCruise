# LesionTOADS volume smoothing-energy fixture

This fixture uses coordinate-distinct float32 class, lesion, and black-hole
memberships on a `5×4×6` volume with `OBJECT`, `LESION`, and `BLACKHOLE`
classes. Per-voxel/class data-term and initial-energy maps are deterministic
binary fractions. A selected-method Java formula harness applies the six-face
smoothing equations at every interior location/class in the source loop order
and records the 11 energy/neighbor results as float32 bit patterns.

The fixture validates streaming composition over fixed inputs. It does not
mutate memberships, apply atlas priors, normalize, or establish lesion
segmentation equivalence. Regenerate with
`scripts/generate_java_lesion_toads_volume_smoothing_reference.sh`.

Selected source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
