# TOADS memory investigation — 2026-10-03

## Trigger and target

The user reported host memory usage of approximately 340 GB against a 40 GB machine limit, and noted that the original Java MaCRUISE run completes at brain resolution below 16 GB. The Python implementation must be brought into that range before resuming long or brain-resolution runs. The 340 GB figure is a user-observed warning; it was not independently measured by this process.

## Findings in the Python implementation

1. `toads_propagate_growing` appended a `ToadsGrowingQueueStep` for every heap event. Each step retained full segmentation, skeleton, availability, ordering, and class-volume arrays. This makes retained memory proportional to the number of growth events times the image-state size.
2. Both growth and thinning queue transitions defensively copied full image arrays for each selected heap item. The thinning neighbor scheduler copied the full availability array after removals as well. Even where old copies were short-lived, this caused large transient memory pressure and allocator churn.
3. The registered loop accumulated all prior exact-iteration and registered-cycle objects, which retain image-sized membership and field arrays.

These paths explain a credible unbounded multiplier; this is a source inspection diagnosis, not proof that any single one caused the full reported peak.

## Changes

- Queue drivers now update their own mutable working arrays in place. Public single-transition helpers retain defensive-copy behavior by default.
- `toads_propagate_growing(record_events=False)` returns no event snapshots and reuses one state. It preserves a trace-recording mode for small source-comparison fixtures.
- Registered loop configuration defaults to no growth-event trace and no prior-iteration history. `iteration_count` reports loop progress when history is disabled. Diagnostic trace scripts explicitly retain history and events.
- The regular CLI reports `iteration_count`, so it does not require keeping all iteration objects alive.

## Validation so far

Command:

```bash
./.venv/bin/pytest -q \
  tests/test_toads_growing_pass.py \
  tests/test_toads_thinning_queue_run.py \
  tests/test_toads_post_registration_branch.py \
  tests/test_toads_registered_first_iteration_cycle.py
```

Result: **9 passed in 17.46s**. The added synthetic 9³ test compares the no-event-history growth path against the Java harness reference for terminal segmentation, skeleton, availability, ordering, and class volumes. Existing queue and registered branch tests also pass.

Follow-up bounded-history and source-reference command:

```bash
./.venv/bin/pytest -q \
  tests/test_toads_registered_first_iteration_cycle.py \
  tests/test_toads_growing_pass.py \
  tests/test_toads_thinning_queue_run.py
```

Result: **6 passed in 23.59s**. On the folded synthetic registered-loop fixture, history-retaining and bounded modes produced the same iteration count, terminal energy ratio, cropped classification, and native classification; bounded mode returned empty history tuples. The base-run restart regression also passed: **2 passed in 74.52s**, exercising the production CLI configuration without iteration/event history.

## Remaining validation before a long run

- Add a direct Java-reference terminal-state comparison for the in-place thinning driver.
- Verify that history-disabled registered-loop results have the same terminal output and iteration count as history-retaining mode on a small synthetic case.
- Profile a controlled brain-resolution synthetic case with an OS-level memory cap and record observed peak RSS. Do not restart the previously interrupted long experiment until this is safe.
- The `<16 GB` target remains unverified for Python until that measured run completes. Synthetic tests establish numerical state equivalence only on their fixtures; they do not establish full MaCRUISE source equivalence or real-data validity.
