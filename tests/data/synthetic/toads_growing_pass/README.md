# TOADS initialized growing pass reference

One deterministic 9×9×9 four-class synthetic label volume contains a folded,
stair-stepped tissue band, a compact second tissue region, and an edge-adjacent
label region. Two skeleton seeds initiate competing fronts. The fixture uses
an asymmetric formula for float32 memberships and ordering values and initial
class volumes `[1, 2, 3, 4]`.

The Java source harness runs the selected `initGrowingLabels` initialization
and `propagateGrowing` loop until the min-heap is empty. It records each popped
priority/payload, each topology result (`null` means stale skeleton), scheduled
neighbors, final segmentation/skeleton/availability/order/volume outputs, and
the residual queue. The harness uses the selected 18/6 LUT and Java
`BinaryHeap4D`. No imaging data are included.

Regenerate with `scripts/generate_java_toads_growing_pass_reference.sh`. This
validates a synthetic growing pass only, not full atlas registration or the
complete TOADS workflow.
