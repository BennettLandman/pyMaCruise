# TOADS outer-loop control fixture

The Java source-formula harness covers the post-membership control from
`AlgorithmMultipleToads`: float32 energy ratio, convergence and iteration-limit
stop decisions, optional growing on the stopping iteration, final
classification update, and relation refresh only when the loop continues.
It includes continued, converged, and iteration-limit cases with both values
of `doGrowing` where relevant.

This validates loop control only. It does not compose the numerical stages or
claim full TOADS segmentation equivalence. Regenerate with
`scripts/generate_java_toads_iteration_control_reference.sh`.

Selected source:
`toads-cruise-scm-2026-09-26/VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java`,
SHA-256 `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`.
