# P18av: TOADS thinning-label and boundary-queue initialization

## Scope

This ports the selected `MultipleToadSegmentation.initThinningLabels` slice:
reset and accumulate float32 class volumes, copy the hard segmentation into
the interior skeleton, mark current-class availability, detect six-face
boundaries, compute source membership-ratio priorities, and seed the P18au
max-heap. It preserves supplied border skeleton and availability state. The
7×7×7 synthetic phantom has four classes; all 125 interior voxels are queued.
No image payload is needed for this formula-level fixture.

The selected source hashes are:

- `MultipleToadSegmentation.java`: `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`
- `BinaryHeap4D.java`: `9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a`

Regenerate the source output with
`scripts/generate_java_toads_thinning_queue_initialization_reference.sh`.
The Java harness calls the retained source `BinaryHeap4D`; the Python package
has no JVM dependency.

## Source endpoint counter finding

The source checks all 26 noncenter neighbors. Every one of the 20 diagonal
neighbors increments `Nngb`; only face neighbors with a different label skip
the increment. Therefore `Nngb` is always at least 20, and the later `Nngb ==
1` endpoint test cannot be reached for an interior voxel. The port retains this
source behavior. The reference queue has no negative-zero endpoint priorities.

## Validation and limits

`tests/test_toads_thinning_queue.py` matches all class-volume float32 bits,
interior skeleton labels, availability flags, per-class boundary counts
`[32, 30, 33, 30]`, and all 125 popped priority-bit/coordinate/class tuples
exactly against Java. Border sentinels and availability values remain
unchanged. The focused queue, heap, and local-transition modules pass **5
tests**.

This validates queue initialization only. Dynamic queue updates after each
thinning decision, whole-volume repeated thinning/growing, topology correction,
atlas deformation, and complete TOADS segmentation remain open.

The full configured suite passes **357 tests in 74.52 seconds** on
implementation commit `ed96f8a3773fe9d5f071b06a51882c0f23dcd23a`; see
[`full-pytest-2026-09-28.md`](full-pytest-2026-09-28.md).
