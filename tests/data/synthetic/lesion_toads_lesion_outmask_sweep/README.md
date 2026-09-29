# LesionTOADS LESION+OUTMASK volume fixture

This synthetic `5×4×6` float32 volume contains ordinary OBJECT classes, one
LESION class with a separate lesion membership/data map, and one OUTMASK class
with a separate outlier membership/data map. It uses mixed direct/prior
relations, nonzero six-face smoothing and atlas-prior terms, asymmetric class
and lesion similarities, zero/unit auxiliary-map values, and zero shape
support. The selected-Java formula harness records every interior class,
lesion, and outlier membership bit after a single source-ordered sweep and
shared normalization.

The fixture validates this selected membership update, not full LesionTOADS
equivalence. The same update with a zero prior coefficient must reproduce its
no-atlas-prior control exactly. Regenerate with
`scripts/generate_java_lesion_toads_lesion_outmask_sweep_reference.sh`.

Selected source:
`toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
