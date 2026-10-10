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

At 2026-10-09 22:07 UTC, preprocessing, crop, initial alignment, registered recrop, and registered-branch checkpoints were complete. The full manifest and lesion subworkflow manifest still reported `running`; the active stage was `registered_loop`. After more than five hours, the stage had produced no checkpoint, per-iteration log, or output file. Peak sampled RSS remained 0.740 GB under the 8 GB watchdog. The prior process handle was later lost during turn interruption; at that point the manifests still showed `running`, but no registered-loop output existed. I restarted the same command with `--resume` at 2026-10-09 23:59 UTC. The watchdog confirmed the resumed process was live at under 0.51 GB RSS in its first minute; it reuses the saved preprocessing, crop, alignment, recrop, and registered-branch checkpoints. This resumed attempt has not yet produced a terminal result. For context, the earlier `plus-retained-center-2026-10-09` synthetic run completed this registered loop at iteration 6, but its tissue labels and ribbon failed the new acceptance checks. The current case has more populated atlas topology classes and its runtime cannot yet be compared at an equal completed-iteration count. The original attempt is a strong practical runtime/observability no-go signal for this run configuration; algorithm correctness and ribbon quality remain undetermined.

This is synthetic workflow execution evidence only. It is not source-equivalence evidence, a completed workflow, or real-MRI validation. No test suite was run specifically for this active-run record. The fixture commit is `ffbdf2b`; this report and its plan/README references are documentation-only.

## Acceptance checks after completion

Run the committed quality audit against the completed terminal workflow manifest:

```sh
./.venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/plus-vandy-atlas-aligned-v3/plus-full-workflow.json \
  --report docs/reports/plus-aligned-folded-ribbon-quality-2026-10-09.json
```

Record and assess RAS grids, discrete labels, watertight connected genus-zero meshes, surface nesting, label support, production triangle closest-point thickness, and level-set subtraction thickness. The existing acceptance bounds are explicit in the audit script. Compare final surface outputs with the P08 Java gates before marking P41z validated. Update this report with the terminal status, elapsed time, results, and the commit that contains the final evidence.
