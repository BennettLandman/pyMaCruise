# LesionTOADS BLACKHOLE+OUTMASK volume fixture

Two synthetic `5×4×6` float32 Java fixtures exercise class orderings with the
BLACKHOLE class before and after OUTMASK. They contain ordinary classes,
separate lesion/black-hole/outlier maps, class and auxiliary data terms,
nonzero class/lesion/black-hole priors, six-face smoothing, mixed relations,
ventricle/GM/brainstem distance maps, and shape gating. The outlier map covers
zero and unit weights. Every interior class and auxiliary membership bit is
recorded after source-ordered updates and shared normalization.

This fixture checks the selected combined membership branch, not full
LesionTOADS equivalence. Both class orders are required because each voxel
update writes memberships and auxiliary maps before processing the next
class. Regenerate with
`scripts/generate_java_lesion_toads_blackhole_outmask_sweep_reference.sh`.

Selected source:
`toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
