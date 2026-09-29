# P18ba — TOADS growing-queue initialization

## Result

Added toads_initialize_growing_queue for the selected initGrowingLabels
behavior. It resets interior availability according to whether the skeleton
is empty, scans voxels and classes in source order, seeds each available class
with a matching six-face skeleton neighbor, computes the float32 priority
ordering + maxMembershipRatio, and returns a min-heap. Existing border state
is preserved; border availability is required to be false because the source
scans border coordinates while reading six neighbors.

Three synthetic 7³ cases match the Java reference:

| Case | Java heap entries | Coverage |
|---|---:|---|
| mixed | 24 | Four adjacent seed classes and multiple candidate classes at shared empty voxels |
| single_seed | 6 | One class-0 seed and its six face candidates |
| no_boundary | 0 | No empty voxels adjacent to an available skeleton label |

For each case, Python matches the Java skeleton, all availability values, and
the complete float32 priority/payload min-heap pop sequence. Membership and
ordering values are asymmetric deterministic synthetic arrays.

## Validation evidence

| Check | Result |
|---|---|
| Java heap-reference cases | 3 cases; 24, 6, and 0 entries; exact arrays and min-heap order |
| Focused pytest | 1 passed in 2.41 s |
| Full pytest | 362 passed in 83.34 s |
| Tested implementation commit | 4be9627c27617db8b6aa94fe3155eda4c71d1718 |
| Full-suite command | .venv/bin/pytest -q from pymacruise/ |
| Environment | pytest 9.1.1, Python 3.13.2 |
| Run date | 2026-09-29 |
| Java harness SHA-256 | f2b5d501a665408e4684af1248f1c625c86cac73aa896cd219bfc34701350ee1 |
| Fixture SHA-256 | 29862a324ff2a616a1cac885ba0c0e093cec1c209f0e5a625da52d91d0c3432a |

Regenerate with
[scripts/generate_java_toads_growing_queue_initialization_reference.sh](../../scripts/generate_java_toads_growing_queue_initialization_reference.sh).
The test and reference fixture are
[tests/test_toads_growing_queue_initialization.py](../../tests/test_toads_growing_queue_initialization.py)
and
[tests/data/synthetic/toads_growing_queue/java_reference.json](../../tests/data/synthetic/toads_growing_queue/java_reference.json).

## Scope limits

This validates initGrowingLabels only. It does not validate a topology-gated
growth event, repeated growth, or integration with the preceding thinning
pass. No actual imaging data was used.
