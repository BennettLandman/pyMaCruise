# TOADS memory investigation — 2026-10-03

The in-place queue and bounded-history implementation is committed locally as `da691f5`; the whole-brain synthetic fixture generator and environment limitations are recorded in `0708e76`. The compact thinning-item-history follow-up is in the current working tree. Commits are local and have not been pushed.

## Trigger and target

The user reported host memory usage of approximately 340 GB against a 40 GB machine limit, and noted that the original Java MaCRUISE run completes at brain resolution below 16 GB. The Python implementation must be brought into that range before resuming long or brain-resolution runs. The 340 GB figure is a user-observed warning; it was not independently measured by this process.

## Findings in the Python implementation

1. `toads_propagate_growing` appended a `ToadsGrowingQueueStep` for every heap event. Each step retained full segmentation, skeleton, availability, ordering, and class-volume arrays. This makes retained memory proportional to the number of growth events times the image-state size.
2. Both growth and thinning queue transitions defensively copied full image arrays for each selected heap item. The thinning neighbor scheduler copied the full availability array after removals as well. Even where old copies were short-lived, this caused large transient memory pressure and allocator churn.
3. The registered loop accumulated all prior exact-iteration and registered-cycle objects, which retain image-sized membership and field arrays.
4. The thinning queue retained one Python heap-item object per popped voxel. These objects are much smaller than full image snapshots, but their count also grows with the number of boundary voxels.

These paths explain a credible unbounded multiplier; this is a source inspection diagnosis, not proof that any single one caused the full reported peak.

## Changes

- Queue drivers now update their own mutable working arrays in place. Public single-transition helpers retain defensive-copy behavior by default.
- `toads_propagate_growing(record_events=False)` returns no event snapshots and reuses one state. It preserves a trace-recording mode for small source-comparison fixtures.
- Registered loop configuration defaults to no growth-event trace and no prior-iteration history. `iteration_count` reports loop progress when history is disabled. Diagnostic trace scripts explicitly retain history and events.
- Production registered runs also omit per-item thinning history. The exact scalar `processed_count` is preserved; standalone kernels continue recording queue items by default for source comparisons.
- The regular CLI reports `iteration_count`, so it does not require keeping all iteration objects alive.
- The synthetic Vandy atlas fixture generator now accepts `--target-shape` and `--output-dir`. A whole-brain-sized `(172, 216, 164)` phantom was generated under ignored `build/vandy-brain-resolution-memory/`; the T1, brain mask, MaSeg, and topology files load as RAS. The compressed generated directory is 52 MB. This avoids committing a second copy of bundled atlas-derived test volumes.

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

After adding compact thinning queue results, the combined focused command passed **7 tests in 91.50s** across `test_toads_thinning_queue_run.py`, `test_toads_default_first_thinning.py`, `test_toads_registered_first_iteration_cycle.py`, and `test_base_workflow_restartable.py`. For each Java queue-run case, the history-free driver preserved terminal segmentation, skeleton, availability, volumes, processed count, and exact remaining queue while returning no popped-item tuple. After strengthening the assertion to compare remaining heap entries, `pytest -q tests/test_toads_thinning_queue_run.py` passed (**1 passed in 0.50s**).

## Whole-brain synthetic profiling run (in progress)

macOS rejected an `RLIMIT_AS` reduction (`ValueError: current limit exceeds maximum limit`) and `/usr/bin/time -l` could not read `kern.clockrate` (`Operation not permitted`). I then implemented a direct macOS `libproc` RSS watchdog. A controlled over-limit test allocated a 200 MB bytearray under a 0.10 GB threshold; the watchdog terminated the child at a sampled 0.209 GB RSS. This confirms that the monitor detects sustained growth and also demonstrates polling overshoot. It is a sampled per-process watchdog, not a hard allocation limit.

The whole-brain-sized synthetic base run is in progress with the watchdog set to 16 GB and a 0.05 second sample interval:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.05 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 build/vandy-brain-resolution-memory/t1.nii.gz \
  --brain-mask build/vandy-brain-resolution-memory/brain-mask.nii.gz \
  --maseg build/vandy-brain-resolution-memory/maseg.nii.gz \
  --atlas-descriptor build/vandy-brain-resolution-memory/atlas.txt \
  --output-dir build/vandy-brain-resolution-memory/profile-run-20261003 \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --toads-iterations 2 \
  --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1
```

At the latest status sample, N4 had completed, `toads-first` was active, and the Python child RSS was 4.19 GB. This is a live sample rather than the run peak. The whole-brain synthetic command has not completed yet; do not claim under-16-GB peak until the watchdog emits its final sampled RSS. `RLIMIT_AS` and `/usr/bin/time` findings explain why the watchdog was added; they are no longer blockers to observation.

The `<16 GB` target remains unverified for Python until the run completes. Synthetic tests establish numerical state equivalence only on their fixtures; they do not establish full MaCRUISE source equivalence or real-data validity.
