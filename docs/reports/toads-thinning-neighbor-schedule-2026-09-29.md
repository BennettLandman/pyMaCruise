# P18aw: TOADS post-removal neighbor scheduling

## Scope

This ports the six-face neighbor-scheduling loop called after a thinning
queue item is removed. The API consumes post-removal segmentation, current
availability, memberships, a max heap, the removed voxel/class, popped
priority, post-removal class volume, volume floor, and limit. It evaluates face
neighbors in source order, calculates the source membership-ratio score (or
uses the popped priority for the source's one-count branch), applies strict
`classVolume > minVolume` and `order > -limit` checks, then marks and queues
accepted candidates.

The Java formula harness uses six deterministic 5×5×5 cases and the actual
selected Java `BinaryHeap4D`. The VandyCRUISE `MultipleToadSegmentation.java`
source hash is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`; the
heap source hash is
`9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a`.
Regenerate the references with
`scripts/generate_java_toads_thinning_neighbor_schedule_reference.sh`.

## Source neighbor-count behavior

The inner loop condition uses the outer face-offset values `i/j/l`, not its
inner offsets `ni/nj/nl`. Because a face offset is always nonzero, the Java
method counts every value in the candidate's full 3×3×3 segmentation window.
The port deliberately preserves this behavior. This count differs from the
separate `initThinningLabels` boundary/endpoint counter documented in P18av.

## Validation and limits

`tests/test_toads_thinning_neighbor_schedule.py` compares exact inserted
priority bits, availability mutations, and full resulting heap pop order. The
cases cover all six directions, source count equal to one, equal priorities,
availability gating, `volume == minVolume`, `order == -limit`, and a priority
just above the strict limit. The focused queue/heap/transition modules pass
**6 tests**.

This validates the neighbor-scheduling primitive only. It does not repeatedly
pop the queue, perform whole-volume thinning/growing, run the topology
correction stages, deform the atlas, or establish final TOADS equivalence.

The full configured suite passes **358 tests in 83.24 seconds** on
implementation commit `f1305fd76ef0a597f90782ac691eb6a53e73a445`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md).
