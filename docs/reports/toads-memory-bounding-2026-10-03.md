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

Four whole-brain-sized (172, 216, 164) base attempts reached the TOADS-first stage but were interrupted before the stage returned. The first attempt had no live periodic output and was stopped when process monitoring became unavailable; attempts 02 and 03 were stopped after about 25 minutes; attempt 04 had live phase logs and was stopped after about an hour in the first registered-loop unity-field/thinning startup. The sampled RSS peak on attempt 04 was **4.299 GB**. N4 completed and is reusable from the manifest, but these incomplete attempts do not establish a full-run peak or completion. Attempts are preserved in separate directories. The RSS watchdog now reports periodically, kills the child process group on interruption or monitor failure, and passed its cleanup test (**1 passed in 0.17s**).

This attempt establishes that the monitored process had remained well below 16 GB in the observed samples; it does not establish the true peak, full-run peak, or completion. The previous watchdog only reported at exit, which made a long TOADS stage hard to observe. `RLIMIT_AS` and `/usr/bin/time` findings explain why the watchdog was added; they are no longer blockers to observation.

### Attempt 05: queue-preflight fix and exact-membership runtime

After moving queue input finiteness checks to the run-level preflight, attempt 05 completed registered unity-field/thinning startup in **523.6 seconds** (TOADS progress time `+488.6s` to `+1012.2s`). The preceding attempt had spent about an hour in this same phase without completing. This confirms that the repeated per-event O(volume) validation was a major startup runtime defect. At the end of startup, sampled RSS was 3.664 GB.

The run then entered `continue_registered_unity_startup_through_first_exact_iteration`. Source inspection showed that its fixed-input membership pass still performs nested Python loops over x, y, z, 11 classes, six neighbors, and other classes, preserving Java's Gauss-Seidel order. On the (172, 216, 164) volume, this is not a useful integration runtime in the current interpreter. I interrupted the run cleanly after roughly two additional minutes in that first exact-membership update; the watchdog reported a peak sampled RSS of **3.903 GB**, well below its 16 GB limit. The manifest keeps the incomplete attempt; only N4 is complete. This run is not a completed whole-brain workflow and does not establish full-run memory or numerical equivalence.

To address this runtime path without changing update ordering, `toads_default_atlas_membership.py` now extracts the source-ordered float32 sweep as a kernel that is Numba-JIT compiled when Numba is installed. The pure-Python implementation remains the fallback for small fixtures. Numba is declared as a runtime dependency, but it is missing from the local environment and package-index access failed earlier, so the compiled kernel has **not** yet been exercised. On commit `139e238`, the fallback passed **13 focused membership and registered-loop tests in 23.65s**:

```bash
./.venv/bin/pytest -q tests/test_toads_default_atlas_dispatch.py tests/test_toads_default_11_membership_continuation.py tests/test_toads_registered_initial_membership.py tests/test_toads_registered_first_iteration_cycle.py
```

These establish unchanged current-interpreter results only, not JIT parity. Next: run these tests under Numba, compare all returned arrays/objective against source harness fixtures, and rerun the whole-brain synthetic workflow with the 16 GB watchdog.

Earlier incomplete attempts established only that sampled per-process RSS remained at or below 4.014 GB before interruption. The completed run below supersedes their whole-workflow runtime and RSS status; the incomplete-attempt history remains useful for explaining how the bottleneck was isolated. Synthetic component tests establish numerical state equivalence only on their fixtures; the full route below does not establish whole-MaCRUISE source equivalence or real-data validity.

### Vectorized initialization and small Vandy-atlas workflow

Source inspection found that thinning initialization visited every interior voxel twice in nested Python loops and checked the 26-neighborhood even for non-boundary voxels. `toads_initialize_thinning_queue` now vectorizes interior skeleton/availability initialization, source-order float32 class-volume accumulation, and six-face boundary detection. It still inserts queue items in Java x/y/z order, preserving equal-priority heap behavior. A literal small-grid Java-order oracle and selected Java harness cases verify skeleton, availability, exact float32 class volumes, boundary counts, and popped queue sequence, including nonpositive/tied membership values.

The focused validation set passed **21 tests in 175.85s** across vectorized initialization, selected TOADS registration/loop, CLI, and base restart behavior. After adding assertions that the base-run log retains live TOADS progress messages, `pytest -q tests/test_base_workflow_restartable.py` passed (**2 passed in 63.98s**). The direct edge-case literal-loop test passed (**1 passed in 0.35s**).

The full-volume runtime bottleneck was traced to repeated validation: `toads_process_next_thinning_item` scanned all memberships and transformed shape values for finiteness on every popped queue item, and neighbor scheduling scanned memberships again for each eligible removal. These O(volume) checks were moved to one-time `toads_run_thinning_queue` preflight validation; public single-item and neighbor APIs still validate their inputs by default. The selected Java queue/startup cases passed after this change (**9 passed in 8.60s**), preserving their terminal arrays, processed counts, and queue references.

The Vandy-atlas (36, 45, 34) synthetic base run completed all six stage checks under the 16 GB RSS watchdog. Peak sampled RSS was **0.598 GB**. TOADS-first and TOADS-second each stopped at the configured two-iteration limit without reaching the energy threshold, as expected for this short integration setting. The manifest marks the run complete but non-converged; all final NIfTI products were verified RAS. A second invocation with `--resume` verified/skipped all six stages. Run artifacts are under ignored `build/vandy-brain-resolution-memory/downsampled-base-run-20261003/`; the runnable command is in `readme.md`.

The Vandy-atlas (36, 45, 34) run exposed the completed path and reported 20.0 seconds for unity-field thinning startup. At brain scale, attempt 04 entered the same startup but ran for about an hour without the completion callback; its peak RSS remained 4.299 GB. This makes repeated O(volume) scans per queue event a source-confirmed runtime defect and a strong explanation for the discrepancy, but the next monitored brain-scale run must verify the effect of moving those scans to one-time preflight. The downsampled completion and incomplete brain-size samples do not establish full Java equivalence or real-data validity.

### Completed brain-resolution synthetic base run

With Numba 0.68.0 installed in the ignored repository `.venv`, the compiled source-ordered membership and statistics paths completed the synthetic base route at `(172, 216, 164)` input resolution. This resolves the earlier runtime blocker for a short, end-to-end synthetic run and confirms the queue-preflight speedup at this scale. It is still a smoke configuration: registration uses one level/iteration, both TOADS passes are capped at two iterations, Dura uses one iteration, and ACE diffusion/TGDM each use one iteration.

The input fixture is generated by `tests/data/synthetic/toads_downsampled_vandy_atlas/generate_fixture.py` with target shape `(172, 216, 164)` into ignored `build/vandy-brain-resolution-memory/`. To repeat from the repository root with a fresh output directory:

```bash
./.venv/bin/python tests/data/synthetic/toads_downsampled_vandy_atlas/generate_fixture.py \
  --target-shape 172 216 164 \
  --output-dir build/vandy-brain-resolution-memory
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.05 --report-interval-seconds 20 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 build/vandy-brain-resolution-memory/t1.nii.gz \
  --brain-mask build/vandy-brain-resolution-memory/brain-mask.nii.gz \
  --maseg build/vandy-brain-resolution-memory/maseg.nii.gz \
  --atlas-descriptor build/vandy-brain-resolution-memory/atlas.txt \
  --output-dir build/vandy-brain-resolution-memory/fresh-run \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --toads-iterations 2 \
  --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --tgdm-iterations 1
```

The first monitored invocation completed with exit code 0. Manifest `run.json` records all six stages as complete, `status=complete`, `converged=false`, `source_equivalence_claim=false`, and the workflow scope as synthetic/native-grid. The TOADS passes each stopped at their two-iteration limit: first pass terminal energy ratio **1.10341454**, second pass **0.570566177**; neither reached the configured energy threshold. Native-grid final products are `inner.nii.gz`, `central.nii.gz`, `outer.nii.gz`, and `SEG_refine.nii.gz`, all on shape `(172,216,164)` with NIfTI axis codes `RAS`. The CRUISE intermediates and four-dimensional GVF also passed their stage check with RAS headers. The manifest records hashes for products and stage logs.

Peak sampled RSS over the complete run was **6.615 GB**, below the 16 GB watchdog threshold; the watchdog did not terminate it. A second invocation with `--resume` reported “verified …; skipping” for N4, TOADS-first, Dura, TOADS-second, CRUISE core, and MaSeg refinement, exited 0, and peaked at **0.123 GB**. These observations are per-process sampled RSS; the watchdog is not a hard OS memory limit and may overshoot between samples.

Observed TOADS phase telemetry: TOADS-first stage completed in **681.7 s** and TOADS-second in **649.0 s**. In each registered update, after the first exact-membership pass, the first atlas registration took **3.2 s**; the growth queues processed **2,052,818 events in 54.4 s** (first pass) and **1,937,385 events in about 35 s** (second pass). This contrasts with earlier whole-brain attempts stalled during repeated per-item whole-volume validation. Whole-run wall time was not recorded by the runner and is left unspecified.

The latest focused checks on the implementation were: JIT exercised in selected Java-reference membership tests (**13 passed**); ordered centroid/statistics references (**26 passed in 6.71 s**); multichannel registration-handoff regressions (**20 passed in 17.86 s**); progress-telemetry plus registered-cycle/base restart regressions (**11 passed in 48.63 s**); and growth queue references after one-time preflight (**3 passed in 0.97 s**). Implementation commits are `139e238`, `8aea116`, `75f5e8b`, and `17fc960`. No full suite has been run after these changes; the older full-suite result does not validate this exact implementation. The next integration work is to match TOADS source initialization grids with Java, establish convergence without smoke caps, then carry that resolved base route into plus workflow. Real-MRI validation remains deferred.

### Follow-up: compile the atlas-support bounds scan

The brain-size run showed that each TOADS start spent 125–131 seconds preparing cropped intensities and initial memberships after atlas arrays loaded. Code inspection found `find_toads_crop_bounds` used `numpy.argwhere` to allocate all foreground coordinates and then mapped them one at a time in Python. That scalar coordinate transform is now an `njit(cache=True, fastmath=False)` kernel with explicit float32 rounding at the same coordinate and matrix-accumulation steps; the Python implementation remains available when Numba is absent.

Focused Java crop-geometry validation after this change:

```bash
./.venv/bin/pytest -q tests/test_toads_crop.py
```

Result: **4 passed in 1.13 s**. The test asserts the compiled dispatcher is exercised when Numba is available, and the pinned Java harness crop bounds remain exact. A supplementary 72³ folded synthetic atlas mask benchmark produced identical bounds from compiled and Python kernels: **0.0013 s compiled versus 0.5324 s Python (419×)** after compilation. This is a microbenchmark, not a brain-resolution rerun; no whole-run speedup is claimed. The benchmark array was generated in memory and was not saved as a test fixture.

This code change postdates the complete base workflow evidence above, which remains tied to `17fc960`. A direct two-iteration TOADS stage was then run with the `e64383c` implementation against the same whole-brain inputs, N4 output, atlas descriptor, and selected registration settings:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.05 --report-interval-seconds 30 -- \
  ./.venv/bin/pymacruise toads-loop \
  --t1-n4 build/vandy-brain-resolution-memory/profile-run-20261003/stages/n4/attempt-01/t1-n4.nii.gz \
  --brain-mask build/vandy-brain-resolution-memory/brain-mask.nii.gz \
  --output build/vandy-brain-resolution-memory/cropjit-probe/classification.nii.gz \
  --dura-output-dir build/vandy-brain-resolution-memory/cropjit-probe/dura-inputs \
  --atlas-descriptor build/vandy-brain-resolution-memory/profile-run-20261003/resources/dura-atlas/atlas.txt \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --maximum-iterations 2
```

It exited 0 after **633.0 s**, stayed under the watchdog with **5.941 GB** peak sampled RSS, and reached the configured iteration limit without reaching the energy threshold. The resulting uint8 classification has shape `(172,216,164)`, axis codes RAS, and SHA-256 `29755467c8835c36f7926a5925332fb39cc9cb421197341021f73c9da50de0c3`; it is **bitwise identical** to the classification from the earlier two-iteration whole-brain run. Crop plus initial-membership setup completed at **111.8 s**, versus **125.9–130.8 s** in prior logs. This is one run per implementation and includes other machine/setup variation; attribute no more than the focused benchmark result to the kernel itself. A requested one-iteration probe was rejected by the CLI, which requires at least two; it wrote no outputs.

The direct probe validates the changed crop code in the whole-brain input geometry, not complete base workflow restartability under the new commit.

### Six-iteration brain-resolution TOADS trace

To measure early convergence behavior at brain resolution without repeating Dura and CRUISE stages, `toads-loop` was run directly with the same N4 T1, mask, packaged Vandy atlas descriptor, and one-level/one-iteration registration settings. Only the TOADS iteration cap changed, from two to six:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.05 --report-interval-seconds 30 -- \
  ./.venv/bin/pymacruise toads-loop \
  --t1-n4 build/vandy-brain-resolution-memory/profile-run-20261003/stages/n4/attempt-01/t1-n4.nii.gz \
  --brain-mask build/vandy-brain-resolution-memory/brain-mask.nii.gz \
  --output build/vandy-brain-resolution-memory/toads-six-iterations-20261003/classification.nii.gz \
  --dura-output-dir build/vandy-brain-resolution-memory/toads-six-iterations-20261003/dura-inputs \
  --atlas-descriptor build/vandy-brain-resolution-memory/profile-run-20261003/resources/dura-atlas/atlas.txt \
  --registration-levels 1 --registration-iterations 1 \
  --registration-max-iterations 1 --maximum-iterations 6
```

The run exited 0 after **1,138.8 s**, with peak sampled RSS **7.199 GB** and no watchdog termination. It completed six iterations but remained nonconverged at the configured cap. The same-input two-iteration run had terminal energy ratio `1.10341454`; the six-iteration run ended at `0.123219952`. Across the six-iteration trace, **1,272,698** classification voxels differed from initialization. The final compressed NIfTI is uint8, `(172,216,164)`, RAS, labels `[1,5,6,14,15,16,17,18,22,24,25]`, SHA-256 `00511035744d02bf1055958422ca1a1d3b2d1fbf125795b4b03830a5a12e04fe`. This provides early synthetic convergence/runtime evidence only; no source whole-loop comparison or clinical interpretation is claimed. The end-of-run ratio is lower than in the two-iteration run but remains well above `1e-4`.

The next integration work is to compare the per-iteration states and objective calculation with the selected Java source at a source-aligned grid, then continue a longer brain-resolution run if the source comparison supports the same trajectory. Memory stayed below the user-reported 16 GB target for this observed process; no real MRI was used.

### Ten-iteration brain-resolution TOADS trace

The same direct `toads-loop` command and inputs were run with `--maximum-iterations 10`, writing to `build/vandy-brain-resolution-memory/toads-ten-iterations-20261003/`. The full command differs from the six-iteration reproduction above only in output directory and iteration cap. It exited 0 after **1,567.7 s**, with **8.053 GB** peak sampled RSS, below the 16 GB watchdog. It completed all ten iterations but remained above the `1e-4` convergence threshold: terminal energy ratio **0.0216144342**, with `iteration_limit_reached=true`. The same-input two-iteration trace ended at `1.10341454`; the six-iteration trace ended at `0.123219952`. These are endpoint measurements from separate deterministic runs, not logged per-iteration ratios.

The final classification is uint8 with shape `(172,216,164)`, RAS axes, labels `[1,5,6,14,15,16,17,18,22,24,25]`, and SHA-256 `1b4f18ad0094c145523f0900ae0777db011be7548cb1e09f5e1380e97e64b615`. This is synthetic execution evidence; it does not establish Java whole-loop equivalence or real-data validity. The observed runtime is about 26 minutes for the TOADS stage under the reduced one-level/one-iteration registration profile. Longer source-default runs may take materially longer. The ratios decrease across the measured endpoints but remain above threshold, so this does not yet validate convergence.

The restartable workflow now records each continuing registered iteration's energy ratio and stop flag in the TOADS stage log. Validation:

```bash
./.venv/bin/pytest -q tests/test_base_workflow_restartable.py
```

Result: **2 passed in 39.20 s**. The test checks that the per-iteration telemetry appears alongside the final run manifest's stop metadata. This focused integration result does not replace a full-suite run.

The full repository pytest suite subsequently passed **1,138 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 793.64 s** on `ad34321`, with peak sampled RSS **13.960 GB** under the 16 GB watcher. The pytest process exited successfully; the watcher reported an error in its final poll/cleanup race, fixed and covered by **2 passed in 0.21 s** on `8190a5d`. The full-suite record is [here](full-pytest-2026-10-03-brain-resolution.md). The selected TOADS Java source-line estimate was re-evaluated after that run and remains 3,464–4,864 of 13,651 lines (25–36%), with zero new source methods counted.

### Compact convergence and registered-loop timing

After adding elapsed-time fields to registered-loop progress messages, the
restartable workflow test passed **2 tests in 20.38 s** on `cdc3f0d`. A direct
two-iteration `toads-loop` on the compact atlas-derived `(36,45,34)` fixture
reported 0.807 s for unity-field/thinning startup, 0.079 s for the first exact
membership iteration, 2.832 s for the first registered cycle, and 4.383 s for
iteration 2. It completed in 10.0 s total, peaked at 0.531 GB sampled RSS, and
ended at energy ratio `0.248976111` (iteration cap, not converged).

The same compact input and selected registration profile were extended to a
maximum of 30 iterations. The loop stopped at iteration 12 after **38.8 s**;
terminal ratio was `7.53300919e-5`, below the configured `1e-4` threshold.
Peak sampled RSS was **0.593 GB**. This is synthetic grid-specific convergence
evidence only.

### Fifteen-iteration brain-resolution continuation

The direct brain-resolution `toads-loop` run used the existing N4 T1, supplied
matching mask, bundled Vandy atlas descriptor, one registration level and one
registration iteration, and a 15-iteration cap. It ran on timing-instrumented
implementation commit `cdc3f0d` under the 16 GB RSS watchdog and completed
without watchdog termination. Total time was **2,103.2 s**; peak sampled RSS
was **7.995 GB**. Phase measurements were 82.420 s for unity-field/thinning
startup, 2.850 s for the first exact-membership iteration, 135.850 s for the
first registered cycle, and 109.842–115.004 s for iterations 3–14. Iteration
15 took 48.388 s because it took the terminal branch.

Per-iteration energy ratios were:

| Iteration | Energy ratio | Iteration | Energy ratio |
|---:|---:|---:|---:|
| 2 | 1.10341454 | 9 | 0.0289222747 |
| 3 | 0.710631251 | 10 | 0.0216144342 |
| 4 | 1.0652374 | 11 | 0.232217193 |
| 5 | 0.205469877 | 12 | 0.131759584 |
| 6 | 0.123219952 | 13 | 0.0726132318 |
| 7 | 0.100688219 | 14 | 0.0362984054 |
| 8 | 0.0455383174 | 15 | 0.0139147937 |

The result is **nonconverged** at the 15-iteration cap. It is a uint8 RAS
volume of shape `(172,216,164)` with labels `[1,5,6,14,15,16,17,18,22,24,25]`
and compressed-file SHA-256
`8673a38edec9ee25b3b5ec75d91333330e48f6084a6ad893043422173f36ecee`. Ratios
are deterministic through iteration 10 relative to the earlier run; iterations
11–15 are new evidence. The nonmonotone trajectory means the compact-grid
convergence result cannot be generalized to brain resolution. The 7.995 GB
observed peak is within the user's 16 GB target, but is a sampled process RSS,
not a hard allocation limit.

### Compact Java/Python registered-loop comparison

The first comparison mistakenly let Python perform five post-refresh rigid
updates while the Java synthetic harness used `alignIter=1`. That transform
comparison was not valid. The Python composition default now matches the
selected Java default of one update, and `toads-loop` and `base-run` expose
`--post-refresh-registration-iterations` so a run can state this control
explicitly.

With matched settings, Java and Python use the same native `(36,45,34)` grid,
internal `(52,60,51)` grid, crop bounds `(2,33),(3,42),(0,30)`, masked N4
input, full-resolution bundled atlas, and one post-refresh registration
update. The initial transform is the same. The post-refresh transform is
`[0.00089594,0.00089594,0.00089594,0.00089594,-2.249104,-6.782193]` in Java
and `[0.00089622,0.00089622,0.00089622,0.00089622,-2.249104,-6.782192]` in
Python (maximum parameter difference below `5e-7`). The transformed
shape-prior arrays differ by at most `7.4e-6` (mean absolute difference
`3.8e-8`). Constructor segmentation, closest-label map, and the initial
membership input to this registration match byte-for-byte; initial-membership
SHA-256 is
`82ec8b65e4c332653ac1011ccd46a436b30b5819f46e439dc52192faaadfc9a0`.

The first captured divergence occurs after the first main-loop thinning,
statistics, and membership update. The Java `computeMemberships()` membership
array and Python's first exact-iteration membership array differ at 23,240 of
55,080 voxels when a voxel is counted if any class membership differs by more
than `1e-6`; maximum absolute channel difference is `0.41345`, and mean
absolute difference over all class values is `0.0006444`. The energy ratios
are `0.60919476` in Java and `0.656713724` in Python. This indicates the
next seam to isolate is first-loop thinning versus the following
centroid/variance/membership recomputation, rather than post-refresh
registration.

After two loop iterations, mapped ordinal-label voxel agreement is
54,011/55,080 (**98.059%**). Per-class Dice values are:

| Class ordinal | Dice |
|---:|---:|
| 1 | 0.995775 |
| 2 | 0.923919 |
| 3 | 0.944072 |
| 4 | 0.825593 |
| 5 | 0.961166 |
| 6 | 0.985915 |
| 7 | 0.836066 |
| 8 | 0.962343 |
| 9 | 0.928571 |
| 10 | 0.930876 |
| 11 | 0.969930 |

Only classes 1 and 6 meet the user's 0.97 minimum Dice. This selected
compact synthetic case is not whole-workflow equivalence, convergence evidence
at brain resolution, or real-data validation. The full brain-resolution
15-iteration trace remains nonconverged at `0.0139147937` after 2,103 seconds
and peaked at 7.995 GB sampled RSS under the 16 GB watchdog.

All exploratory harness files, converted raw arrays, logs, and images are
inside ignored `build/java-reference/` or `build/downsampled-vandy-base-run/`;
no source repository files outside `pymacruise` were changed. No plus-workflow
or real-MRI validation was attempted in this cycle.

Next: compare the second rigid-registration coefficients, transformed-prior input, and line-search trajectory against the selected Java implementation. Once that seam is aligned, rerun the compact and brain-resolution convergence traces.
