# Current-code 64×80×64 synthetic MaCRUISEplus workflow — 2026-10-08

## Scope and provenance

This reruns the naturally converged, restartable raw-input plus workflow on the committed 64×80×64 folded synthetic T1, FLAIR, user brain mask, and MaSeg inputs, using the current implementation fingerprint `b200db7a08b9bec962249dab5de65aa0be11f650324a8e8699b65e42256a748f`. The previous code generation was preserved under the same output root; the current run is isolated under `code-revisions/b200db7a08b9bec9/`.

This is synthetic integration and resource evidence. It does not claim full MATLAB/JIST composite equivalence or validate a human MRI. The plus-resolution comparison still has a known sensitivity to input-grid discretization; see [the grid surface scaling report](plus-grid-surface-scaling-2026-10-06.md).

## Reproduction

```sh
env PYTHONNOUSERSITE=1 PYTHONPATH=src \
  MPLCONFIGDIR=build/wheel-plus-natural/mpl-cache \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_scaled_grid/t1.nii.gz \
  --flair tests/data/synthetic/plus_scaled_grid/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_scaled_grid/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_scaled_grid/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt \
  --output-dir build/plus-scaled-grid-current/natural-run-64 \
  --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
  --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
  --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
  --loop-thinning-limit 1000 --finite-zero-variance --resume
```

The existing complete manifest belonged to an older implementation fingerprint. With `--resume`, the runner created a new implementation-specific stage root and ran the current code while preserving that prior generation. It exited 0, naturally converged, and the RSS watchdog did not terminate it; peak sampled RSS was **0.824 GB** under the 16 GB limit.

## Stage and output evidence

The top-level manifest `build/plus-scaled-grid-current/natural-run-64/plus-full-workflow.json` reports `status=complete` and `terminal=true`. All nine stages completed: preprocessing, crop, initial alignment, registered recrop, LesionTOADS, Dura products, native-grid transfer, MaSeg cortical separation, and the plus surface workflow. LesionTOADS stopped by convergence at **iteration 6** with energy ratio **`7.722568398094154e-07`**, below `1e-4`; the iteration cap was not reached.

The complete set of nested stage manifests contains **72 unique hashed output records**: 62 compressed NIfTIs, 3 VTK meshes, and 7 summaries/provenance files. An independent audit verified every path/hash record across the new stage root, loaded every output NIfTI, and confirmed finite voxels with RAS orientation. It parsed each VTK surface and confirmed finite points, valid triangles, and in-range point indices. The top-level plus workflow records 17 final product links; the downstream surface stage includes inner/central/outer level sets, three RAS-mm surfaces, and thickness volume/metrics.

## Restart and tests

The same command and settings were rerun with `--resume` and a stricter 4 GB watchdog. The workflow reported that it verified and resumed the raw-input plus LesionTOADS workflow, exited 0, and peaked at **0.236 GB RSS**.

The full repository suite on code-equivalent commit `5a412b7` passed **1,243 tests, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 1,743.86s**, with **14.961 GB** peak sampled RSS; see [the suite report](full-pytest-2026-10-08-base-surface-thickness.md). The full suite was not rerun for this workflow execution. This run advances synthetic plus workflow integration only; plus source-composite comparison, the cross-resolution positional tolerance issue, and real-MRI validation remain open.
