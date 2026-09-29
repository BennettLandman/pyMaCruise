# TOADS growing queue one-event reference

Three deterministic 5×5×5 cases compare one `propagateGrowing` queue event
with the selected VandyCRUISE Java source and its `BinaryHeap4D` and 18/6
critical-point LUT: a topology-accepted event that schedules five exposed
face-neighbors, a topology-rejected event, and a stale queue item whose
skeleton center is already occupied. Memberships use the established integer
formula divided by 30 in float32. Volumes begin at `[10, 20, 30, 40]`; the
initial min-heap contains priorities `0.25` (event) and `0.5` (remaining item).

The reference records the popped item, topology result, complete segmentation,
skeleton and availability arrays, center ordering value, float32 class-volume
bits, and complete remaining Java heap pop order. This synthetic source-slice
test does not validate repeated growth or end-to-end TOADS. There is no subject
data.

Regenerate with
`scripts/generate_java_toads_growing_queue_step_reference.sh`.
