# Current-source synthetic MaCRUISEplus workflow — 2026-10-08

## Scope

This is a fresh, restartable, 64×80×64 synthetic plus run using the current VandyCRUISE pyramid implementation (`069e9d3`, implementation signature `abe8c0ecb3a5127e7540d515d345ca425f458d35508f00eeab4e9b722611f457`). It validates current-code workflow integration and the registered LesionTOADS convergence path. It does not claim whole-workflow MATLAB/JIST equivalence or real-MRI validation.

## Reproduction

```sh
env PYTHONNOUSERSITE=1 PYTHONPATH=src \
  MPLCONFIGDIR=build/plus-pyramid-source-current/mpl-cache \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_scaled_grid/t1.nii.gz \
  --flair tests/data/synthetic/plus_scaled_grid/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_scaled_grid/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_scaled_grid/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt \
  --output-dir build/plus-pyramid-source-current/natural-run \
  --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
  --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
  --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
  --loop-thinning-limit 1000 --finite-zero-variance
```

The run exited 0 under the 16 GB watchdog; peak sampled RSS was **0.791 GB**. All nine top-level stages completed. LesionTOADS stopped by convergence at **iteration 6** with energy ratio **`7.722568398094154e-07`**, below the configured `1e-4` threshold; the 99-pass cap was not reached.

An independent recursive audit of every `outputs` and `final_outputs` reference across the nested JSON manifests verified **75 unique files**, comprising **62 compressed NIfTIs**, **3 VTK surfaces**, and **10 ancillary files** (summaries, manifests, and the ANTs transform). All NIfTIs contain finite values and are RAS-oriented. The workflow manifest lists 17 final product links.

## Restart and interpretation

The same command and settings with `--resume` exited 0 and verified/resumed the raw-input plus LesionTOADS workflow. Peak sampled RSS was **0.108 GB**. This confirms current-source synthetic integration, convergence, output integrity, and restart behavior. It does not establish end-to-end source equivalence. To reconcile counts, the prior report listed 72 principal records (62 NIfTIs, 3 meshes, and 7 selected ancillary outputs). Applying the recursive all-reference audit to that preserved implementation generation verifies 74 unique files; the current generation has 75, adding the Dura output checkpoint manifest. All 62 NIfTIs validate in RAS in both generations. The source-level comparison for the changed pyramid method is documented separately in [the VandyCRUISE actual-method report](toads-registration-pyramid-actual-source-2026-10-08.md).
