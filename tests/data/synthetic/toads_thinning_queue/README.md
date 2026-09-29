# TOADS thinning-label/queue initialization phantom

This 7×7×7, four-class array phantom is generated from deterministic integer
formulas; it contains no subject images. The Java reference initializes the
skeleton, current-class float32 volumes, availability flags, boundary counts,
and the max-heap in selected `MultipleToadSegmentation.initThinningLabels`
source order. It compares all 125 queued interior items by float32 priority
bits and coordinate/class payload, plus all output arrays. Border sentinel and
availability values are verified unchanged.

The selected method's source counter increments for every diagonal
neighbor, so its `Nngb == 1` endpoint branch cannot be reached for an interior
3×3×3 neighborhood. The port preserves this source behavior; the phantom's
125 priorities contain no negative-zero endpoint entries. This is a recorded
source quirk, not a corrected interpretation.

Regenerate the reference with
`scripts/generate_java_toads_thinning_queue_initialization_reference.sh`.
The Java class reference uses `BinaryHeap4D.java` from the test-only source
copy, not a runtime JVM.
