# LesionTOADS three-auxiliary volume fixture

Two synthetic `5×4×6` float32 Java formula fixtures exercise one each of
LESION, BLACKHOLE, and OUTMASK, with the selected `lesionClass` being BLACKHOLE
or LESION according to which special class appears last in atlas class order.
They cover separate lesion/black-hole terms and maps, class/lesion/black-hole
similarity priors, both distance-modulation branches, OUTMASK-specific
neighbor energy and outlier blending, direct/prior relations, zero shape
support, and source normalization of the shared lesion map at both special
classes. All interior class, lesion, black-hole, and outlier membership bits
are stored.

This validates the selected multi-class formula branch, not full LesionTOADS
equivalence. Both atlas orders are necessary because class updates are
in-place. Regenerate with
`scripts/generate_java_lesion_toads_three_auxiliary_sweep_reference.sh`.

Selected source:
`toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
