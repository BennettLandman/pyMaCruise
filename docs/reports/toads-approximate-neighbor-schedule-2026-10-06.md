# Approximate TOADS neighbor scheduling — 2026-10-06

## Scope

Added `schedule_approximate_thinning_neighbors` for the six-face rescheduling loop at the end of the selected VandyCRUISE `MultipleToadSegmentation.propagateThinning` event. It carries sparse class slots into the source `maxApproxMembershipRatio`, preserves nested face-neighbor traversal, float32 priority arithmetic, the class-volume gate, availability mutation, and the selected source's literal 26-neighbor count condition.

## Evidence

- Synthetic case: the committed folded 9×9×9 three-class/two-slot fixture from `tests/test_toads_approximate_thinning_step.py`.
- Selected Java source SHA-256: `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
- Standalone Java excerpt harness SHA-256: `483a3de5f045d093cd5a84978fd0c425415d2847b15daed0d6e8462a13cace9e`.
- Pinned event: priority bits `-1082130432`, coordinate `(1,1,1)`, class `0`.
- Java and Python insert the same two neighbors in source order, each with priority bits `-1073741824`: `(1,2,1,0)` and `(2,1,1,0)`.
- The fixture opens the class-volume gate with a synthetic zero minimum and generous priority limit. It is event-level source evidence, not realistic minimum-volume behavior or approximate-path convergence.
- Focused command: `./.venv/bin/python -m pytest -q --tb=short tests/test_toads_approximate_thinning_step.py` — **2 passed in 0.55s**.
- Full suite on implementation commit `3f37a7e`: **1,204 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,214.15s**, peak sampled RSS **14.419 GB** under the 24 GB watchdog; see [full-suite report](full-pytest-2026-10-06-approx-neighbor-schedule.md).

## Remaining approximate-path work

The complete queue lifecycle still needs its production minimum-volume threshold behavior, repeated event processing, approximate centroid/deviation updates, membership recalculation, registration control, growing, stopping, and output transfer composed and validated. This row does not claim a complete `useApprox` constructor run.
