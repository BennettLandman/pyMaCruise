# TOADS binary heap source reference

The JSON fixture records four deterministic operation sequences against the
selected Java `BinaryHeap4D`: min and max ordering with duplicate priorities
and forced capacity growth, mode change with entries already queued, and reset
followed by reuse. Values are compared by float32 bit pattern and each popped
`(x, y, z, label)` payload is exact. The cases include signed zero and
coordinates/labels that narrow through Java signed-short conversion.

The test-only source copy retains its July 2004 version and Pierre-Louis Bazin
author comments at
`tests/java_reference/upstream/edu/jhmi/rad/medic/structures/BinaryHeap4D.java`.
Its SHA-256 is
`9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a`.
Regenerate the fixture with
`scripts/generate_java_toads_binary_heap_reference.sh`.

This validates the queue primitive only; use in TOADS thinning/growing and
complete segmentation remains separate work.
