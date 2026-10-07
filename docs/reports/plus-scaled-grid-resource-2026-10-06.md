# MaCRUISEplus 64×80×64 resource investigation — 2026-10-06

## Purpose

P41w measures whether the naturally converging raw-input plus workflow scales
from the committed 36×45×34 folded synthetic phantom to the larger committed
64×80×64 synthetic phantom. This is synthetic workflow/resource testing; no
human MRI was used.

Fixture generator and provenance are committed under
`tests/data/synthetic/plus_scaled_grid/`. Inputs are T1, FLAIR, user brain mask,
and MaSeg, all compressed NIfTI and RAS. The run uses the selected mini
LesionTOADS atlas descriptor.

## Runs and observations

The first installed-wheel `plus-full-run` attempt used the 24 GB RSS watchdog.
It reached the registered LesionTOADS branch, but a `proc_pidinfo` lookup
failure caused the watchdog to stop the child. The output manifests remained
`running`, so that attempt is inconclusive and was not counted as a workflow
failure.

The watchdog now tolerates a transient process lookup failure when the child is
still alive, retries briefly, and fails closed after repeated failures or when
the prior sampled peak approaches the configured ceiling. Focused regression
command:

```text
./.venv/bin/python -m pytest -q --tb=short tests/test_rss_watchdog.py
```

Result: **3 passed, 0 failed, 0 skipped in 1.94s** on commit
`7e06817cb23f7835c82892ad5432916448f83e97`. The regression covers a child
exiting between the poll and RSS read, and a one-time lookup failure while the
child remains alive.

The complete installed-wheel workflow was then resumed with a 32 GB ceiling,
leaving 8 GB below the user's stated 40 GB machine limit. It reached the same
registered LesionTOADS branch. The watchdog sampled **25.825 GB peak RSS**;
the child then exited with status **-9** while `terminated_for_rss_limit=false`.
The run did not produce the registered-branch output manifest, terminal plus
outputs, or a successful resume. Cause of the `-9` exit is unconfirmed; it may
reflect an external process or execution-environment limit. It is not evidence
of an algorithm exception or a successful workflow. The memory figure is the
watchdog's sampled child RSS and may miss brief peaks.

## Instrumented retry and registration hot path

The registered-branch checkpoint now records atomic substage markers in its
manifest. A source-checkout retry under the 24 GB ceiling showed this sequence:
initial state hydration, statistics refresh, distance-map refresh, and
GM-enabled membership all completed; `shape_registration:running` remained
active. After roughly ten minutes without leaving that stage, the run was
interrupted deliberately through the watchdog so the identified slow stage
could be improved. It sampled **0.329 GB peak RSS** and did not reach growth,
relation refresh, or output materialization. This did not reproduce the earlier
25.825 GB peak; the location of that peak remains unknown.

The selected registration path called both TOADS and LesionTOADS shape
samplers at every sampled voxel, and each sampler rescanned an entire atlas
volume for finite values even though the registration-volume caller had
already validated those maps. That repeated full-map validation is now skipped
only through an internal prevalidated-input path. The standard-13 membership
pass also computed the same centroid-similarity matrices at every voxel; these
are now computed once before the Java-ordered scan, as invariant iteration
state. Commit `c437b238a624effdd50a8556a9445ec9d45c5e76` passed the source-order
Java standard-13 membership harness and focused shape/checkpoint regressions:
**21 passed, 0 failed, 0 skipped in 27.19s**, peak sampled RSS **0.460 GB** under
an 8 GB watchdog. This confirms fixture-level arithmetic but does not yet
measure the 64×80×64 optimized registration stage.

Reproduction command (run from `build/wheel-plus-natural/`, using the isolated
installed package at `site/`):

```text
env PYTHONNOUSERSITE=1 PYTHONPATH=/Users/landmaba/dev/mipav/pymacruise/build/wheel-plus-natural/site MPLCONFIGDIR=/Users/landmaba/dev/mipav/pymacruise/build/wheel-plus-natural/mpl-cache /Users/landmaba/dev/mipav/pymacruise/.venv/bin/python /Users/landmaba/dev/mipav/pymacruise/scripts/run_with_rss_watchdog.py --limit-gb 32 --interval-seconds 1 --report-interval-seconds 60 -- /Users/landmaba/dev/mipav/pymacruise/.venv/bin/python -m pymacruise.cli plus-full-run --t1 /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/plus_scaled_grid/t1.nii.gz --flair /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/plus_scaled_grid/flair.nii.gz --brain-mask /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/plus_scaled_grid/brain-mask.nii.gz --maseg /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/plus_scaled_grid/maseg.nii.gz --atlas-descriptor /Users/landmaba/dev/mipav/pymacruise/tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt --output-dir natural-run-64 --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 --loop-thinning-limit 1000 --finite-zero-variance --resume
```

## Interpretation and next work

The 64×80×64 run is **not validated**. The partial checkpoint is retained under
the ignored `build/wheel-plus-natural/natural-run-64/` directory. Do not count
the nine-stage plus workflow as complete on this grid or infer brain-resolution
capacity from it. P41x is measuring the optimized registration and later
branch stages under bounded RSS while preserving selected source order and
numerical outputs. The acceptance target is a complete,
naturally converged run, finite RAS outputs, verified hashes and resume, plus a
sampled RSS peak below a conservative 24 GB ceiling.

## Follow-up on 2026-10-06: compact growth traces and restart

A source-checkout run at `build/plus-scaled-grid-current/natural-run-64/`
completed the initial first-iteration checkpoint, shape registration, and all
13 shape-bound calculations. It entered `growth:running`, sampled a **20.221 GB
peak RSS** under a 24 GB watchdog, then exited `-9` with
`terminated_for_rss_limit=false`. The exact external kill trigger remains
unconfirmed.

Code inspection found that `propagate_lesion_toads_growing` retained a full
copy of segmentation, skeleton, availability, ordering, and class-volume
arrays in every per-pop event. This made trace storage scale as queue length
times volume size. Commit `2e0c39f` adds compact event records for integrated
volume runs while keeping full event snapshots as the default for source
comparison callers. The compact mode preserves popped queue items, topology
decisions, and inserted neighbors, and final image arrays are unchanged.
Focused checks on commit `2e0c39f`:

```text
./.venv/bin/python -m pytest -q tests/test_lesion_growing_pass.py tests/test_plus_lesion_registered_branch.py tests/test_rss_watchdog.py
```

Result: **13 passed in 142.10s**, peak sampled RSS **0.439 GB** under an 8 GB
watchdog. A focused checkpoint/resume rerun after allowing a missing downstream
registered-loop manifest to be created on parent resume passed **2 tests in
26.67s**, peak sampled RSS **0.441 GB**. In the 64-grid workflow, resumed first
iteration and registered-branch checkpoints were hash-verified and reused. The
registered branch then completed growth, relation refresh, and output
materialization at **0.523 GB peak RSS** under the 24 GB ceiling.

The top-level workflow next failed because its not-yet-created registered-loop
child manifest was rejected when resuming the parent. That restart edge is now
fixed: a missing downstream child starts from its completed upstream manifest.
A second `--resume` invocation is currently running the registered-loop stage
with compact growth traces. Its latest sampled peak is **0.574 GB**; it has not
yet produced the loop manifest or final Dura/CRUISE outputs. P41w remains **not
validated** until natural convergence, terminal RAS outputs, hashes, resume,
and the remaining plus stages all complete on the 64 grid.

The full raw-input command writes only inside the ignored `build/` tree. The
latest retry uses the existing `natural-run-64` directory and `--resume`; see
`run-resume2.log` for RSS samples. Keep the 24 GB watchdog ceiling for this
test. Do not count the nine-stage plus workflow as complete on the 64 grid or
infer brain-resolution capacity before terminal outputs are validated.
