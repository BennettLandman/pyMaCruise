# LesionTOADS OUTMASK atlas-prior volume fixture

This synthetic `5×4×6` float32 volume extends the OBJECT+OUTMASK sweep fixture
with nonzero, asymmetric class-shape priors and off-diagonal class similarities.
It uses mixed direct/prior relations, an outlier map with zero/interior/unit
weights, and a zero-shape OUTMASK voxel. The selected-Java formula harness
records every interior class and outlier membership bit after smoothing, atlas
prior accumulation, inversion, OUTMASK blending, and normalization.

The fixture validates this selected OBJECT+OUTMASK branch, not full
LesionTOADS equivalence. Setting `prior_coefficient=0` must reproduce the
P39zw smoothing-only result. Regenerate with
`scripts/generate_java_lesion_toads_outmask_prior_sweep_reference.sh`.

Selected source:
`toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
