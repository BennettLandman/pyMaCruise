# TOADS relation-state refresh fixture

This committed fixture records a selected-Java composition of the
`AlgorithmMultipleToads` refresh calls on one asymmetric 9×8×7 folded
categorical phantom. It carries a nonuniform float32 class-membership tensor
and prepopulated signed-byte boundary-distance state through
`computeRelations`, `computeBoundaryDistances`, and `computeClosestLabels`, in
that order. The JSON stores every output map/table/slot/byte and the Java loop
counts. There is no medical image data.

Regenerate with
`scripts/generate_java_toads_relation_state_refresh_reference.sh`. The source
file is `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`,
SHA-256 `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
This validates only the composed selected methods, not TOADS segmentation
evolution or approximate-membership mode.
