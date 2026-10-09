# Corrected-code synthetic MaCRUISEplus rerun

## Outcome

The folded 64×80×64 synthetic plus workflow completed all nine stages on
implementation commit `b785bff` (`Match retained TOADS registration center`).
The registered LesionTOADS loop naturally stopped at iteration 6 with energy
ratio `7.722568398094154e-07`; the iteration limit was not reached. The top
manifest records `status=complete` and `terminal=true`.

An independent audit verified all 92 unique file paths with SHA-256 records
referenced by the nested manifests (including inputs and outputs). All 67
referenced NIfTIs were finite and RAS-oriented. The three algorithmic VTK
surface products declare RAS coordinates in millimetres. A same-settings
`--resume` verified and resumed the completed workflow successfully under a
4 GB watchdog, peaking at 0.237 GB. The initial run peaked at 0.800 GB under a
16 GB limit and took approximately 51 minutes from output-directory creation
to final manifest completion; most of the elapsed time was the registered
LesionTOADS loop, which currently provides no per-iteration checkpoint while
it runs.

This confirms synthetic implementation and restart behavior after the shared
TOADS retained-center correction. It does not establish whole-workflow
Java/MATLAB equivalence or real-MRI validation. Per-iteration progress and
runtime optimization for the registered LesionTOADS loop remain open.

## Reproduction

```sh
env PYTHONNOUSERSITE=1 PYTHONPATH=src \
  MPLCONFIGDIR=build/plus-retained-center-2026-10-09/mpl-cache \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_scaled_grid/t1.nii.gz \
  --flair tests/data/synthetic/plus_scaled_grid/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_scaled_grid/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_scaled_grid/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt \
  --output-dir build/plus-retained-center-2026-10-09 \
  --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
  --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
  --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
  --loop-thinning-limit 1000 --finite-zero-variance
```

The resume command repeats the same arguments and appends `--resume`, using a
4 GB watchdog. It reports `Verified and resumed raw-input plus LesionTOADS
workflow` and exits successfully.

The full repository suite on the same implementation commit passed **1,503
tests, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 1,764.57
s**, with **13.787 GB** peak sampled RSS; see [full-suite evidence](full-pytest-2026-10-09-retained-center.md).
