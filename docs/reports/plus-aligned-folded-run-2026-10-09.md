# Atlas-aligned folded synthetic plus run

## Purpose

Exercise the raw-input `plus-full-run` workflow on a compact atlas-aligned synthetic case with distinct CSF, white-matter, and gray-matter MaSeg labels, then apply explicit ribbon-quality checks to the resulting production surfaces. This report distinguishes workflow progress from a completed validation result.

## Inputs and command

The committed fixture is `tests/data/synthetic/plus_vandy_atlas_aligned/`, a shared RAS grid of 65×81×62 voxels. Its MaSeg image contains discrete values 0, 5, 44, and 100; its topology image contains 13 discrete classes from the bundled LesionTOADS atlas. T1, FLAIR, brain mask, MaSeg, and atlas descriptor are all NIfTI or text assets within the repository. The fixture generator is `scripts/generate_plus_vandy_atlas_aligned.py`.

The run was started with:

```sh
env PYTHONNOUSERSITE=1 \
  PYTHONPATH=src \
  MPLCONFIGDIR=build/plus-vandy-atlas-aligned-v3/mpl-cache \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
    --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
    ./.venv/bin/python -m pymacruise.cli plus-full-run \
      --t1 tests/data/synthetic/plus_vandy_atlas_aligned/t1.nii.gz \
      --flair tests/data/synthetic/plus_vandy_atlas_aligned/flair.nii.gz \
      --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned/brain-mask.nii.gz \
      --maseg tests/data/synthetic/plus_vandy_atlas_aligned/maseg.nii.gz \
      --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned/cruise-atlas-12obj-lesiontoads2012.txt \
      --output-dir build/plus-vandy-atlas-aligned-v3 \
      --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
      --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
      --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
      --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
      --loop-thinning-limit 1000 --finite-zero-variance
```

## Current state

The first attempt was still in `registered_loop` after more than five hours and produced no loop checkpoint; peak sampled RSS was 0.740 GB under the 8 GB watchdog. A resume at 2026-10-09 23:59 UTC likewise did not advance to a loop output, so it was interrupted with saved upstream checkpoints intact. Commit `1751b45` adds a hash-associated JSON progress log written after every completed iteration, including iteration number, energy ratio, stop reason, stage order, and elapsed time. Focused regression validation on that implementation passed **13 tests, 0 failed, 0 skipped in 117.00s**, peak RSS **0.430 GB** under a 6 GB watchdog.

A new full-workflow attempt under this changed implementation signature recalculated the registered branch instead of reusing the older branch checkpoint. Its live operation history showed the delay was in the source-selected `shape_optimizer`. I stopped that duplicate work and started a direct one-iteration continuation from the existing complete branch checkpoint at `build/plus-vandy-atlas-aligned-v3/lesion_workflow/registered-branch/lesion-registered-branch.json`, writing to `build/plus-vandy-atlas-aligned-v3/loop-chunks/iteration-0002`. The progress log was created but had no completed-iteration event before interruption; RSS stayed below 0.7 GB. This gives a concrete runtime bottleneck and improves future observability, but the newly labeled case still has no terminal loop result or ribbon audit. The run cannot yet be compared at equal completed-iteration counts with the earlier `plus-retained-center-2026-10-09` case, which stopped at iteration 6 but failed tissue-label and ribbon-quality checks.

Commit `4ccf6e6` removes repeated temporary class arrays from the source-ordered coefficient reducer. The LesionTOADS registration Java-reference group passes **15/15 in 2.19s**, and the registered-branch integration group passes **12/12 in 101.38s**; both ran under a 6 GB watchdog with peak RSS below 0.43 GB. A follow-up bounded continuation still did not finish its first pass after several minutes. I also evaluated a Numba compiled sampling kernel; it changed exact Java coefficient totals and optimizer trial traces, so that experiment was reverted and is not part of the implementation. Throughput remains an open blocker.

The recorded runtime is a practical throughput blocker for this synthetic case, not a source-equivalence or ribbon-quality verdict. Current workflow evidence is synthetic execution only; it is not source-equivalence evidence, a completed workflow, or real-MRI validation. Fixture commit: `ffbdf2b`; iteration progress instrumentation and focused tests: `1751b45`.

Follow-up on 2026-10-10: a bounded continuation from the saved iteration-1 registered branch completed the source-selected registered-loop pass and wrote a complete, hash-checked checkpoint at `build/plus-vandy-atlas-aligned-v3/loop-chunks/iteration-0003/lesion-registered-loop.json`. It advanced to iteration 2, retained the source stage order, reported energy ratio `0.026117630302906036`, and did not stop; elapsed time was `509.09s`, with peak sampled RSS `0.470 GB`. This is the first completed registered LesionTOADS iteration on the properly labeled plus fixture, but it is an intermediate checkpoint and does not yet permit the terminal surface audit. The run used implementation commit `d61c21e`; 29 focused lesion shape and TOADS sampling tests passed in `2.20s`, including exact equality between direct and reused shared-coordinate sampling. Commit `bcbd7db` then factored the registered-shape affine coordinate once per voxel across classes; its selected initial-thinning, repeated-iteration, and recrop integration tests passed `3/3` in `12.11s`. A continuation from iteration 2 using that optimization is now running, with a separate one-pass checkpoint directory.
That continuation completed at implementation commit `bcbd7db`: iteration 3, energy ratio `0.17062699794769287`, stop=false, elapsed `534.25s`, peak sampled RSS `0.480 GB`. The checkpoint is `build/plus-vandy-atlas-aligned-v3/loop-chunks/iteration-0004/lesion-registered-loop.json`. The iteration-2 to iteration-3 segment was not faster than iteration 1 to 2, so no end-to-end speedup is claimed. The next run will record timings around transformed-prior sampling, the continuing segmentation/statistics/membership pass, and registration/branching before further optimization.
Subsequent phase-timed checkpoints on the same fixture show practical progress. Iteration 4 completed on `90e12dc` in `513.99s` (prior sampling `13.22s`, continuing segmentation/statistics/membership `135.04s`, registration/branching `365.73s`; energy ratio `0.1678145826`, stop=false). A source-ordered slab-vectorized coefficient path, committed as `3992264`, passed the exact Java volume fixtures and reduced the next pass (iteration 5) to `360.48s` (12.91s / 136.79s / 210.78s; ratio `0.1523705572`, stop=false). A vectorized transformed-support-bounds path, committed as `362cb76`, passed selected Java handoff bounds and reduced iteration 6 to `309.87s` (13.03s / 139.50s / 157.32s; ratio `0.1507368386`, stop=false). Peak sampled RSS stayed at or below `0.475 GB` under the 8 GB watchdog. The current checkpoint is `build/plus-vandy-atlas-aligned-v3/loop-chunks/iteration-0007/lesion-registered-loop.json`. This is a 40% total pass-time reduction from the 514s starting measurement; registration/branching time is down about 57%. The continuing segmentation/statistics/membership phase now dominates, and subphase timing is the next instrumented target.

Across completed iterations 2 through 6, all five checkpoints are intermediate, synthetic, and hash-verified; none is a terminal result. Energy ratio rose overall from `0.0261176303` at iteration 2 to `0.1507368386` at iteration 6, after peaking at `0.1706269979` at iteration 3; the selected source convergence stop has not occurred. No ribbon-quality or full workflow-equivalence verdict is yet available.

These are synthetic execution and source-reference test results, not full workflow equivalence or real-MRI validation. The ribbon audit remains pending until the selected loop reaches a terminal result and the complete workflow manifest is available.

The fixture generator, full-run invocation, progress-log schema, and eventual audit command remain the reproducibility record. The audit is intentionally pending until a complete terminal manifest exists.

## Acceptance checks after completion

Run the committed quality audit against the completed terminal workflow manifest:

```sh
./.venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/plus-vandy-atlas-aligned-v3/plus-full-workflow.json \
  --report docs/reports/plus-aligned-folded-ribbon-quality-2026-10-09.json
```

Record and assess RAS grids, discrete labels, watertight connected genus-zero meshes, surface nesting, label support, production triangle closest-point thickness, and level-set subtraction thickness. The existing acceptance bounds are explicit in the audit script. Compare final surface outputs with the P08 Java gates before marking P41z validated. Update this report with the terminal status, elapsed time, results, and the commit that contains the final evidence.
