# LesionTOADS one-event growth reference

The 5³ synthetic reference isolates one selected
LesionToadSegmentation.propagateGrowing queue pop. The Java formula harness
uses the checked 18/6 CriticalPointLUT, source Numerics.max, and source
BinaryHeap4D; the topology predicate reuses the separately validated
selected multi-label method body. It covers an accepted lesion event in sum
and max modes, an accepted non-lesion event, a topology rejection, and a
stale item whose skeleton location was already filled.

For accepted cases the fixture captures exact segmentation/skeleton and
availability fields, ordering and class-volume float32 bits, inserted
neighbors, and the complete resulting heap pop order. Rejection re-enables
the event label's availability without changing volumes or inserting
neighbors; stale events skip the topology check and all state updates. The
reference does not instantiate the complete LesionTOADS class and does not
establish full queue exhaustion or segmentation equivalence.

Regenerate the committed fixture with
scripts/generate_java_lesion_growing_event_reference.sh.
