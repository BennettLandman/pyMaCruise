# P18ay — repeated TOADS thinning queue loop

## Result

Added `toads_run_thinning_queue`, which repeatedly invokes the source-matched
P18ax event and follows VandyCRUISE's outer condition literally:
`isNotEmpty() && (val > -limit)`, where `val` is initialized to float32 zero
and thereafter holds the previous popped priority. Thus the first item at or
below `-limit` is still popped and processed when the prior priority is above
the threshold. It returns the final image state, every processed heap item,
and the residual max heap.

The Java reference has four synthetic cases:

| Case | Processed items | Purpose |
|---|---:|---|
| `drain_with_rescheduling` | 8 | Repeated event state changes, neighbor insertion, stale item, and queue exhaustion |
| `cross_threshold_at_second_item` | 2 | The second item is exactly `-limit`; it is processed before the loop stops |
| `top_priority_equals_threshold` | 1 | The first item is exactly `-limit`; it is processed from the initial `val=0` state |
| `zero_limit_processes_nothing` | 0 | Initial strict comparison fails, leaving the heap untouched |

For every processed item, the fixture records the exact priority bits and
payload plus post-event segmentation and skeleton snapshots. It also records
the final availability state, float32 class volumes, and residual heap order.
The queue drain matches the Java outputs exactly.

The harness uses the selected Java `CriticalPointLUT` and `BinaryHeap4D`, and
reuses the event implementation in
`ToadsThinningQueueStepReferenceHarness.java`. VandyCRUISE
`MultipleToadSegmentation.java` SHA-256:
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.

## Validation evidence

| Check | Result |
|---|---|
| Java loop cases | 4 cases; all processed item tuples, event snapshots, final arrays, float32 volumes, and residual heap entries match exactly |
| Focused queue test subset | 4 passed in 0.57 s |
| Full pytest | 360 passed in 84.50 s |
| Tested implementation commit | `858651f71030628bd4adb5ae8db6f950c41b0a60` |
| Full-suite command | `.venv/bin/pytest -q` from `pymacruise/` |
| Environment | pytest 9.1.1, Python 3.13.2 |
| Run date | 2026-09-29 |
| Java harness SHA-256 | `2b55b46bd879b5737b143ff9521911d8ffb41ee765b2d8d576f16a79645b5dd4` |
| Fixture SHA-256 | `7a5a2557d2ba67f24b5cd552eec8a74d40228f7a61d8883384ffd22f9cd7fdf2` |

Sources and regeneration:
[`tests/test_toads_thinning_queue_run.py`](../../tests/test_toads_thinning_queue_run.py),
[`tests/data/synthetic/toads_thinning_queue_run/java_reference.json`](../../tests/data/synthetic/toads_thinning_queue_run/java_reference.json),
and
[`scripts/generate_java_toads_thinning_queue_run_reference.sh`](../../scripts/generate_java_toads_thinning_queue_run_reference.sh).

## Scope limits

This composes the outer loop with one already validated event. It seeds its
synthetic starting state from the P18ax fixture; it does not validate
`initThinningLabels` plus the whole loop as a single operation. It does not
implement whole TOADS segmentation, deformable atlas registration, other
topology-correction passes, or Dura orchestration. No actual imaging data was
used. P18 and P34 remain open.
