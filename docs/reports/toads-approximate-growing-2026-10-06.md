# Sparse approximate TOADS growing pass — 2026-10-06

## Scope

Added `initialize_toads_approximate_growing_queue` and
`toads_approximate_propagate_growing` for the selected
`useApprox` branches in VandyCRUISE `MultipleToadSegmentation.initGrowingLabels`
and `propagateGrowing`. The implementation keeps the sparse candidate/membership
slots and computes `maxApproxMembershipRatio` in float32 source order; it does
not expand memberships to a dense all-class volume. It composes the selected
min-heap, 18/6 topology gate, accepted/rejected voxel updates, sparse class
volume accumulation, ordering state, and six-face neighbor scheduling.

## Source comparison

- Selected VandyCRUISE source: `MultipleToadSegmentation.java`, SHA-256
  `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
- Method slices: `initGrowingLabels` lines 465–500 and `propagateGrowing`
  lines 1783–1856; this validates only the sparse membership arms and their
  queue composition.
- Java harness:
  `tests/java_reference/ToadsApproximateGrowingPassReferenceHarness.java`,
  SHA-256 `fe33597307cd64570f61d95bb0fad58f5cf9357055f728d13441d6e9f4b866a5`.
- Deterministic folded 9×9×9, three-label, two-slot fixture. The complete
  Java/Python event lists match for **480 popped events**, including each
  priority bit pattern, coordinate, class, topology decision, and inserted
  neighbor. There are 293 accepted and 29 rejected topology decisions; the
  other events are stale/non-applicable queue events. Terminal segmentation,
  skeleton, availability, ordering float32 bits, and class-volume float32 bits
  match; the queue drains.
- Focused command:
  `./.venv/bin/python -m pytest -q --tb=short tests/test_toads_approximate_growing.py tests/test_toads_growing_pass.py tests/test_toads_growing_queue_initialization.py tests/test_toads_growing_queue_step.py`
  — **6 passed in 1.71s**.
- Full suite on `9c7a800bbb7b33b774a3a75fdf6e8fb74e48a350`: **1,209 passed,
  0 failed, 0 skipped, 11 expected RAS warnings in 1,309.68s**, peak sampled
  RSS **13.363 GB** under the 24 GB watchdog; see [full-suite report](full-pytest-2026-10-06-approx-growing.md).

This is source-derived synthetic component evidence. It is not the complete
TOADS constructor or iteration loop, and it does not establish Java/MATLAB
workflow equivalence, convergence, or real-data behavior.
