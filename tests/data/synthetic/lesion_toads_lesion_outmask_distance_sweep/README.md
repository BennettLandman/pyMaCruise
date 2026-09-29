# LesionTOADS LESION+OUTMASK distance-volume fixture

This synthetic `5×4×6` float32 fixture extends the selected combined
OBJECT+LESION+OUTMASK membership sweep with signed-byte ventricle, GM, and
brainstem distance maps. Ventricle distances test zero-distance `spread` and
precedence; a `-1` ventricle sentinel exercises GM fallback; brainstem
distances modulate independently, including zero. The maps coexist with
asymmetric atlas priors, mixed direct/prior relations, and varying outlier
weights. The Java harness stores every interior class, lesion, and outlier
membership bit.

The fixture validates selected distance routing in the composed volume
formula, not full LesionTOADS equivalence. Disabling all distance maxima must
reproduce the P39zy combined result exactly. Regenerate with
`scripts/generate_java_lesion_toads_lesion_outmask_distance_sweep_reference.sh`.

Selected source:
`toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
