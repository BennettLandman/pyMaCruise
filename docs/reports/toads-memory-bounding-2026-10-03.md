# TOADS memory investigation — 2026-10-03

The in-place queue and bounded-history implementation is committed locally as `da691f5`; the whole-brain synthetic fixture generator is in `0708e76`; the compact thinning-item-history follow-up is in `b897a45`; and the RSS watchdog plus the active synthetic profile command are in `3c0b132`. Commits are local and have not been pushed.

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

## Whole-brain profiling and completed downsampled run

macOS rejected an `RLIMIT_AS` reduction (`ValueError: current limit exceeds maximum limit`) and `/usr/bin/time -l` could not read `kern.clockrate` (`Operation not permitted`). I then implemented a direct macOS `libproc` RSS watchdog. A controlled over-limit test allocated a 200 MB bytearray under a 0.10 GB threshold; the watchdog terminated the child at a sampled 0.209 GB RSS. This confirms that the monitor detects sustained growth and also demonstrates polling overshoot. It is a sampled per-process watchdog, not a hard allocation limit.

The first whole-brain-sized synthetic base run was started with the watchdog set to 16 GB and a 0.05 second sample interval:

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

Three whole-brain-sized (172, 216, 164) base attempts reached the TOADS-first stage but were interrupted before the stage returned. The first attempt had no live periodic output and was stopped when process monitoring became unavailable; attempts 02 and 03 were explicitly stopped after about 25 minutes without a completed TOADS stage. The final sampled RSS peak across the monitored attempts was **4.014 GB**. N4 completed and is reusable from the manifest, but this is not a full-run peak or completion claim. Attempts are preserved in separate directories. The RSS watchdog now reports periodically, kills the child process group on interruption or monitor failure, and passed its cleanup test (**1 passed in 0.17s**).

This attempt establishes that the monitored process had remained well below 16 GB in the observed samples; it does not establish the true peak, full-run peak, or completion. The previous watchdog only reported at exit, which made a long TOADS stage hard to observe. `RLIMIT_AS` and `/usr/bin/time` findings explain why the watchdog was added; they are no longer blockers to observation.

The `<16 GB` target remains unverified for Python until a whole-brain run completes and its watchdog summary is recorded. Incomplete attempts established only that sampled per-process RSS remained at or below 4.014 GB before interruption; they are not full-run evidence. Synthetic tests establish numerical state equivalence only on their fixtures; they do not establish full MaCRUISE source equivalence or real-data validity.

### Vectorized initialization and small Vandy-atlas workflow

Source inspection found that thinning initialization visited every interior voxel twice in nested Python loops and checked the 26-neighborhood even for non-boundary voxels. `toads_initialize_thinning_queue` now vectorizes interior skeleton/availability initialization, source-order float32 class-volume accumulation, and six-face boundary detection. It still inserts queue items in Java x/y/z order, preserving equal-priority heap behavior. A literal small-grid Java-order oracle and selected Java harness cases verify skeleton, availability, exact float32 class volumes, boundary counts, and popped queue sequence, including nonpositive/tied membership values.

The focused validation set passed **21 tests in 175.85s** across vectorized initialization, selected TOADS registration/loop, CLI, and base restart behavior. After adding assertions that the base-run log retains live TOADS progress messages, `pytest -q tests/test_base_workflow_restartable.py` passed (**2 passed in 63.98s**). The direct edge-case literal-loop test passed (**1 passed in 0.35s**).

The Vandy-atlas (36, 45, 34) synthetic base run completed all six stage checks under the 16 GB RSS watchdog. Peak sampled RSS was **0.598 GB**. TOADS-first and TOADS-second each stopped at the configured two-iteration limit without reaching the energy threshold, as expected for this short integration setting. The manifest marks the run complete but non-converged; all final NIfTI products were verified RAS. A second invocation with `--resume` verified/skipped all six stages. Run artifacts are under ignored `build/vandy-brain-resolution-memory/downsampled-base-run-20261003/`; the runnable command is in `readme.md`.

The full brain-resolution run remains open. A fresh monitored attempt with live stage logs is needed to determine whether the observed long runtime is due to queue initialization, registered TOADS iterations, or their scale dependence. Neither the downsampled completion nor the incomplete brain-size RSS samples establish full Java equivalence or real-data validity.
