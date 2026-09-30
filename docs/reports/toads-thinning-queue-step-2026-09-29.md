# P18ax — one queue-driven TOADS thinning event

## Result

Implemented `toads_process_next_thinning_item` to compose one selected
`propagateThinning` event: pop the max-heap item, apply the source volume and
skeleton guards, choose the best membership class, find a direct or
relation-linked face target, apply the 18/6 homeomorphic LUT, remove the point
from the skeleton, update the selected class volume and availability, and
schedule six-face candidates using the popped priority.

Seven deterministic 5×5×5 cases compare exact segmentation and skeleton arrays,
availability arrays, float32 class volumes, and subsequent heap pop tuples
against the Java reference. They cover direct topology acceptance and
rejection, relation-linked replacement, an unconnected best class, exact
volume-floor rejection, a stale skeleton item, and the source's behavior when
the segmentation label differs from the queue label but the skeleton agrees.
The last case confirms that the popped class `k`, rather than the current
segmentation value, selects the membership and volume update.

The harness runs the selected source `CriticalPointLUT` and `BinaryHeap4D`
classes. The production VandyCRUISE source is
`MultipleToadSegmentation.java`, SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`; the LUT
SHA-256 is
`50f0deaf2328741755cad5cdafa91c0debd5f6a09f5712d52f3945e35e257f9f`; the heap
SHA-256 is
`9d9d016a6d87bd9aa07e77b75429651ee66d4f67422515d5bd8ab3dde6a9ac7a`. The
source method is from April 2006 and names Pierre-Louis Bazin as author.

## Validation evidence

| Check | Result |
|---|---|
| Java LUT/heap reference cases | 7 cases; exact discrete arrays, class volumes, and subsequent heap contents/order |
| Focused pytest | 1 passed in 0.73 s |
| Full pytest | 359 passed in 73.99 s |
| Tested implementation commit | `6b2c74adadf05097a6e21662a6bf7ae07c0c2330` |
| Full-suite command | `.venv/bin/pytest -q` from `pymacruise/` |
| Environment | pytest 9.1.1, Python 3.13.2 |
| Run date | 2026-09-29 |
| Fixture JSON SHA-256 | `c36cebac2325e00de409084b0454da4122ae5d4e9e2af328ec73df1dd3b6888b` |
| Java harness SHA-256 | `e791377912152a1a64be7f4aa89ef467490fced215f41f65f9871574cbb52ffd` |

The test and reference fixture are
[`tests/test_toads_thinning_queue_step.py`](../../tests/test_toads_thinning_queue_step.py)
and
[`tests/data/synthetic/toads_thinning_queue_step/java_reference.json`](../../tests/data/synthetic/toads_thinning_queue_step/java_reference.json).
Regenerate the reference with
[`scripts/generate_java_toads_thinning_queue_step_reference.sh`](../../scripts/generate_java_toads_thinning_queue_step_reference.sh).

## Scope limits

This is one queue event on deliberately small synthetic grids, not repeated
global queue evolution, convergence, deformable TOADS, or a complete
segmentation. It uses the exact-membership (`useApprox == false`) selection
path. No subject imaging data was used. P18 and its P34 dependency remain open.
