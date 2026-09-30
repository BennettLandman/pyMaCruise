# TOADS relation-list fixture

Three synthetic categorical volumes compare the selected
`buildRelationList` source method with Python: uniform labels, a planar two-class
interface, and a folded six-class configuration. References capture every
compact interior relation ID, the constructor-initialized class labels left at
the border, every relation membership bit, and the class-level union rows.
No medical images are used.

Regenerate with `scripts/generate_java_toads_relation_list_reference.sh`.
The test harness transcribes the selected method from
`VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`,
SHA-256 `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
The result validates relation-list construction only, not the subsequent
relation-map diffusion or full TOADS segmentation.
