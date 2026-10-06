# Approximate TOADS queue pass — 2026-10-06

## Scope

Added `run_approximate_thinning_pass`, composing sparse `useApprox` thinning initialization, repeated selected queue-event transitions, six-face neighbor rescheduling, source queue-priority stopping, and class minimum-volume gates. The normal class floor follows `max(minimum_volume_floor, volume_ratio * initial_class_volume)`; explicit class minima are available to make small event fixtures exercise removals. Full array copies and event retention are disabled by default inside the loop, so memory scales with image state rather than processed-event count.

## Evidence

- Synthetic input: deterministic folded 9×9×9 three-class label case with two sparse membership slots; slot 0 equals the current label, so this fixture validates queue lifecycle/removal order without exercising label reassignment/topology branches.
- Selected VandyCRUISE `MultipleToadSegmentation.java` SHA-256: `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
- Java source-derived harness SHA-256: `e15f0898cffb762632ce162b6a2bb1230deff8a5dea9018afac0556aa5536ec6`.
- With explicit zero minimum volumes and generous `limit=1000`, Java and Python process the same **343 events**, drain the queue, match the SHA-256 of the ordered popped priority/coordinate/class trace, and match all three terminal class-volume float32 bit patterns. Both leave no interior skeleton voxels and preserve the label map.
- With the selected default 100-voxel minimum floor, this compact fixture's classes remain below threshold: initial boundary items are popped without removal, matching the source guard's intended behavior.
- Focused command: `./.venv/bin/python -m pytest -q --tb=short tests/test_toads_approximate_thinning_pass.py tests/test_toads_approximate_thinning_step.py tests/test_toads_approximate_thinning_startup.py` — **6 passed in 0.65s**.
- Full suite on implementation commit `2e16f4b`: **1,206 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,212.71s**, peak sampled RSS **14.427 GB** under the 24 GB watchdog; see [full-suite report](full-pytest-2026-10-06-approx-thinning-pass.md).
- Runtime follow-up commit `9f44477` avoids full membership/shape rescans per event; focused queue group: **7 passed in 0.66s**. Full suite on that tested commit: **1,207 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,217.29s**, peak sampled RSS **14.615 GB**; see [validation-boundary full-suite report](full-pytest-2026-10-06-approx-thinning-validation-boundary.md).

## Limits

This is source-derived synthetic queue-stage evidence, not a Java execution of the full constructor or complete `useApprox` TOADS. The fixture does not exercise label reassignment/topology rejection, realistic above-floor brain volumes, approximate centroid and deviation refresh, membership recomputation, rigid registration, growing, convergence, or workflow outputs. These remain separate integration work.
