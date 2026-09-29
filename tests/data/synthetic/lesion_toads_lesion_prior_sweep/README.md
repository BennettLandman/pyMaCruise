# LesionTOADS OBJECT + LESION prior sweep fixture

This synthetic `5×4×6` float32 volume extends the P39zq OBJECT+LESION sweep
with nonzero selected class/lesion shape priors, asymmetric off-diagonal class
similarities, a lesion-similarity vector, mixed direct/prior relations, and a
zero-shape LESION voxel. The selected-Java formula harness follows the
in-place x/y/z update, including the special lesion-class contribution to
ordinary class priors and the lesion-membership prior. It records class and
auxiliary lesion membership float32 bits at every interior voxel.

Distance modulation, black-hole/OUTMASK classes, registration, topology, and
full segmentation equivalence remain outside scope. Regenerate with
`scripts/generate_java_lesion_toads_lesion_prior_sweep_reference.sh`.

Selected source: `toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
