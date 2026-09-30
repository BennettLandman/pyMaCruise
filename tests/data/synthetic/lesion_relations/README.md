# LesionTOADS relation-list and diffusion fixture

Six synthetic label volumes compare Python with selected method excerpts from
`LesionToadSegmentation.java`: an asymmetric planar interface with zero
iterations, a wide planar interface with one and two source limits, an
asymmetric folded multi-class field at limits one and two, and a uniform
field. The wide planar pair exercises propagation over more than one diffusion
front. The JSON records input labels, initial relation map, final
relation map, relation table, relation count, and Java loop count. Cases use
non-cubic dimensions to expose axis-order mistakes. No medical images are used.

The Java excerpt harness implements the selected `buildRelationList` and
`computeRelations` loops and can be regenerated with
`scripts/generate_java_lesion_relations_reference.sh`. The source file is
`VandyCRUISE/src/edu/jhmi/rad/medic/methods/LesionToadSegmentation.java`,
SHA-256 `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`.
`buildRelationList` visits six faces in nested `i/j/l` order; the diffusion
selection order is X-, X+, Y-, Y+, Z-, Z+. Diffusion is synchronous per pass
because Java writes a temporary map and then copies all interior values.

This validates these selected relation methods against the Java excerpt
harness only. It does not validate topology checks, lesion segmentation
updates, or complete LesionTOADS execution.
