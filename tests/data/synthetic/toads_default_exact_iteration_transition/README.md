# Default exact TOADS post-membership transition

This synthetic fixture validates the bounded source-order transition after an
exact default-atlas membership pass. The membership result is represented by
its memberships and objective; tests feed it into the existing topology-aware
growing and relation-refresh implementations. Energies exercise a continuing
iteration and both final-iteration `doGrowing` settings.

The selected VandyCRUISE order in `AlgorithmMultipleToads.runAlgorithm()` is
`computeMemberships()`; optionally `atlas.registerShapes()`; compute the
float32 energy ratio and stop decision; call `propagateGrowing()` when
`(!stop) || doGrowing`; then either call `updateClassificationFromLabels()`
when stopping, or call `computeRelations()`, `computeBoundaryDistances()`, and
`computeClosestLabels()` in that order when continuing.

Atlas registration is held fixed and lies outside this fixture. The arrays
are generated deterministic synthetic maps (the branch tests also use the
separately Java-referenced `toads_growing_pass` phantom); no patient images are
used. This fixture does not establish full TOADS equivalence or validate atlas
registration. Regenerate its Java control reference with
`tests/java_reference/generate_toads_default_exact_iteration_transition_reference.sh`.

Source files and SHA-256 values:

- `VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java`:
  `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`
- `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`:
  `63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884`
