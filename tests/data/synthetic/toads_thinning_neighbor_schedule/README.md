# TOADS post-removal neighbor-scheduling fixture

The six deterministic 5×5×5 cases provide an already updated segmentation,
availability grid, removed label, prior queue entries, source popped priority,
post-removal class volume, minimum volume, and limit. Memberships use the
integer recipe in the Java harness. The harness transcribes the selected
`propagateThinning` six-face scheduling loop and uses the actual test-only
`BinaryHeap4D` implementation.

Cases exercise all six face directions, the source inner-count `== 1` path,
a duplicate-priority heap tie, unavailable-neighbor gating, exact volume
equality, exact `order == -limit`, just-above limit insertion, and final heap
pop order. It preserves a source indexing quirk: the inner count condition
uses the outer face offset, so every face candidate counts all 27 segmentation
samples in its 3×3×3 neighborhood. This is intentionally source-faithful and
kept distinct from the `initThinningLabels` neighbor-count quirk.

Regenerate the JSON with
`scripts/generate_java_toads_thinning_neighbor_schedule_reference.sh`. This is
the post-removal queue update primitive only; it does not run repeated
thinning or TOADS evolution.
