# Repeated TOADS thinning queue-run reference

The fixture records the selected VandyCRUISE outer `propagateThinning` loop
over synthetic state from the committed P18ax queue-step reference. The Java
harness uses the selected source `CriticalPointLUT`, `BinaryHeap4D`, and the
P18ax event harness. Regenerate it with
`scripts/generate_java_toads_thinning_queue_run_reference.sh`.

Cases cover full queue drain with dynamic face-neighbor scheduling, crossing
the negative limit on the second item, popping a top item whose priority is
exactly `-limit`, and zero limit (which processes nothing). They retain the
full sequence of popped priorities and payloads plus segmentation/skeleton
snapshots after every event, final availability and float32 class volumes, and
the remaining heap order.

The two threshold-only cases suppress neighbor scheduling so they isolate the
outer-loop test. The drain case uses `limit=100`, matching the event reference
harness's scheduler cutoff. This validates the repeated event driver and its
source stop rule; the initial state is fixture-backed rather than constructed
by `initThinningLabels` here. It is still not complete TOADS segmentation.
