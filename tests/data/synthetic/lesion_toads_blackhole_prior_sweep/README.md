# LesionTOADS BLACKHOLE prior sweep fixture

The committed selected-Java formula harness extends the `5×4×6` OBJECT plus
BLACKHOLE fixture with asymmetric class priors, off-diagonal class
similarities, lesion/black-hole similarities, mixed relations, and
`priorCoefficient=0.3`. It checks the selected BLACKHOLE class, lesion, and
black-hole prior sums along with source-order smoothing, shape gating,
inversion, and normalization. All interior membership outputs are stored as
float32 bits.

Distance modulation and full segmentation remain excluded. Regenerate with
`scripts/generate_java_lesion_toads_blackhole_prior_sweep_reference.sh`.

Selected source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
