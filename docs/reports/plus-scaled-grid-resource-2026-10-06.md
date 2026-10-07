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

## Follow-up — compact growth traces and 64-grid completion

### Memory cause and correction

A source-checkout run entered registered LesionTOADS growth at 64×80×64,
sampled a **20.221 GB peak RSS** under a 24 GB watchdog, then exited `-9`
without the watchdog terminating it. Inspection found that
`propagate_lesion_toads_growing` retained full copies of segmentation,
skeleton, availability, ordering, and class-volume arrays for every queue
pop. The trace therefore grew with queue length times image size.

Commit `2e0c39f` added a compact event record for integrated volume runs,
retaining each popped item, topology decision, and inserted neighbors while
discarding prior full-volume states. Full snapshots remain the default for
source-comparison callers. The Java-reference folded-growth fixture confirms
identical event traces and final arrays in both modes. The checkpoint resume
path was also repaired so a parent workflow can reuse existing child
manifests and start a downstream stage whose manifest does not exist yet;
commit `1de2109` contains that change.

Focused validation on `2e0c39f`:

```text
./.venv/bin/python -m pytest -q tests/test_lesion_growing_pass.py tests/test_plus_lesion_registered_branch.py tests/test_rss_watchdog.py
```

Result: **13 passed in 142.10s**, peak sampled RSS **0.439 GB** under an 8 GB
watchdog. The focused registered-loop restart check passed **2 tests in
26.67s**, peak sampled RSS **0.441 GB**. No full repository suite was run for
these changes.

### 64×80×64 result

The committed folded synthetic inputs completed all nine raw-input plus
stages and naturally converged at LesionTOADS iteration 6. The energy ratio
was **7.722568398094154e-07**, below the selected `1e-4` tolerance. The
registered branch reached growth, relation refresh, and output materialization
at **0.523 GB peak RSS** after reusing validated first-iteration files. The
registered-loop process peaked at **0.582 GB**. The final downstream Dura,
native transfer, MaSeg, and surface stages completed at **0.389 GB**. The
maximum sampled RSS across the post-fix processes was **0.582 GB**, down from
the preceding 20.221 GB failed attempt.

The final top-level manifest is `complete` and `terminal=true`. All nine stage
manifests completed. Validation checked **71 manifest output hashes**, **61
finite RAS `.nii.gz` images**, and **3 VTK surfaces** declaring RAS coordinates
in millimeters. The final surface thickness summary reports mean **0.0624097
mm** over **3,724 vertices**; volume thickness is available with seven
crossing voxels. A same-settings invocation of
`./scripts/run_plus_scaled_grid_64_synthetic.sh --resume` verified and reused
the completed outputs.

The compact 36×45×34 natural-convergence fixture stops at iteration 5 with
energy ratio `8.241410796472337e-06` and reports mean surface thickness
`0.11254296722799782 mm` over 998 vertices. The 64-grid summary is lower by
about **0.05013 mm** and uses a different voxel grid and mesh tessellation;
these means are not a pointwise or symmetric nearest-surface comparison and
do not establish the agreed 0.01 mm positional acceptance. Composite
Java/MATLAB equivalence and real-MRI validation remain separate.

Summed process time across the four 64-grid invocations was about **2 h 2 min**
(13 min 51 s before the memory fix; 10 min 17 s completing the registered
branch; 1 h 36 min 21 s completing the registered loop; and 1 min 33 s
finishing downstream stages). This is the sum of active command intervals,
not one uninterrupted run; gaps between retries are excluded. The reproducible
command is `./scripts/run_plus_scaled_grid_64_synthetic.sh`; add `--resume` to
reuse completed checkpoints. Outputs stay under the ignored `build/` tree.
No human MRI was used.
