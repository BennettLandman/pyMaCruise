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

The fixture generator, full-run invocation, progress-log schema, and eventual audit command remain the reproducibility record. The audit is intentionally pending until a complete terminal manifest exists.

## Acceptance checks after completion

Run the committed quality audit against the completed terminal workflow manifest:

```sh
./.venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/plus-vandy-atlas-aligned-v3/plus-full-workflow.json \
  --report docs/reports/plus-aligned-folded-ribbon-quality-2026-10-09.json
```

Record and assess RAS grids, discrete labels, watertight connected genus-zero meshes, surface nesting, label support, production triangle closest-point thickness, and level-set subtraction thickness. The existing acceptance bounds are explicit in the audit script. Compare final surface outputs with the P08 Java gates before marking P41z validated. Update this report with the terminal status, elapsed time, results, and the commit that contains the final evidence.
