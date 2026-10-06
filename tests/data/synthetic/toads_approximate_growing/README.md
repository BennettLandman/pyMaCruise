# Sparse approximate TOADS growing phantom

This deterministic 9×9×9 fixture uses three labels with concentric and
folded-cortex-like boundaries. A sparse two-slot candidate/membership array
and a thinned seed skeleton drive the `useApprox` branches of selected
VandyCRUISE `MultipleToadSegmentation.initGrowingLabels` and
`propagateGrowing` (source SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`). The
committed JSON contains the Java harness queue-event trace, including every
inserted min-heap item and topology accept/reject result, plus terminal
segmentation, skeleton, availability, ordering, and float32 volume bits.

Generate the source reference with
`scripts/generate_java_toads_approximate_growing_reference.sh`. The harness
uses the committed selected 18/6 critical-point LUT and Java `BinaryHeap4D`
reference classes. This fixture validates the approximate growing queue
transition on synthetic arrays only; it is not full TOADS or whole-workflow
equivalence.
