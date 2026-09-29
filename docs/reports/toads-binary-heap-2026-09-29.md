# P18au: tie-faithful TOADS four-dimensional heap

## Scope

`src/pymacruise/toads_heap.py` ports the selected VandyCRUISE `BinaryHeap4D`
min/max heap operations for priorities and `(x, y, z, label)` payloads. It
preserves Java float32 priority conversion, signed-short payload narrowing,
strict comparisons during insertion, left-child preference on equal children
during removal, no stable FIFO guarantee, reset, and mode changes without
rebuilding an existing heap. The queue is intended for later TOADS boundary
thinning/growing integration; it does not yet perform that integration.

The original Java class is retained unchanged under
`tests/java_reference/upstream/edu/jhmi/rad/medic/structures/`. Its original
header identifies July 2004 and Pierre-Louis Bazin. SHA-256:
`9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a`.
The Java harness and fixture generator are
`tests/java_reference/ToadsBinaryHeap4DReferenceHarness.java` and
`scripts/generate_java_toads_binary_heap_reference.sh`.

## Validation

`tests/test_toads_heap.py` replays the committed Java scenarios: min/max
ordering with duplicate priorities and capacity growth, changing heap mode
while items remain queued, and resetting/reusing the heap. It compares each
popped float32 bit pattern and coordinate/class payload exactly. The scenarios
include signed zero and Java signed-short narrowing. Two focused tests pass.

The full configured suite passes **355 tests in 73.96 seconds** on
implementation commit `10475c8c2c4f1977c5a9e70b3e9b1bcd649975c3`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md). This is primitive
ordering evidence only, not global TOADS queue initialization or evolution.
