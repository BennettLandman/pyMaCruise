# One-event TOADS thinning reference

This synthetic 5×5×5 fixture records seven one-item queue events executed by
`ToadsThinningQueueStepReferenceHarness.java`. It uses the selected
VandyCRUISE `CriticalPointLUT` and `BinaryHeap4D` sources and the source-ordered
logic from `MultipleToadSegmentation.propagateThinning` (April 2006; author
Pierre-Louis Bazin). Regenerate with
`scripts/generate_java_toads_thinning_queue_step_reference.sh`.

The cases exercise direct topology acceptance and rejection, a relation-linked
replacement, an unconnected best class, exact volume-floor rejection, a stale
skeleton queue item, and the source rule that checks the skeleton and class
volume using the popped queue label even when the current segmentation label
differs. Each case compares the resulting segmentation, skeleton, availability
state, float32 class volumes, and remaining Java heap pop order.

The queue step is a composed source slice, not a full thinning loop. The
reference has two entries and processes exactly one. It does not validate
convergence, repeated scheduling, full-volume TOADS evolution, or a clinical
segmentation. The available face-neighbor bits in each case are seeded
explicitly; all non-center memberships use `[0.1, 0.8, 0.05, 0.05]`, with the
center membership recorded in the JSON input block. Relation and shape-prior
arrays are zero except in `related_allowed`, which enables relation `[2, 1]`
and gives the `(-x)` face neighbor a class-1 shape prior of `0.7`.
