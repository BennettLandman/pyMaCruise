# Fresh synthetic MaCRUISEplus run after TGDM correction — 2026-10-07

## Run and result

A fresh raw-input plus run completed all nine ordered stages on the committed 64×80×64 folded RAS synthetic fixture. It used the T1, shifted-grid FLAIR, matching brain mask, MaSeg labels, and bundled synthetic atlas under `tests/data/synthetic/plus_scaled_grid/` and `tests/data/synthetic/plus_lesion_atlas_init/`. The run did not use a human scan.

```sh
PYTHONNOUSERSITE=1 PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_scaled_grid/t1.nii.gz \
  --flair tests/data/synthetic/plus_scaled_grid/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_scaled_grid/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_scaled_grid/maseg.nii.gz \
  --atlas-descriptor tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt \
  --output-dir build/plus-after-tgdm-2026-10-07/natural-run-64 \
  --registration-iterations 2 --shrink-factors 1 --smoothing-sigmas 0 \
  --registration-sampling-rate 1 --n4-shrink-factor 1 --n4-iterations 2 \
  --n4-spline-mesh 1 1 1 --crop-padding 1 --init-align-iterations 0 \
  --pyramid-levels 2 --lesion-iterations 1 --first-thinning-limit 1000 \
  --loop-thinning-limit 1000 --finite-zero-variance
```

The top-level manifest is `complete` and `terminal=true`. Preprocessing, crop, initial alignment, registered recrop, LesionTOADS, Dura output materialization, native-grid transfer, MaSeg cortical separation, and plus surface/thickness all passed. LesionTOADS converged naturally at iteration **6**, with energy ratio `7.722568398094154e-07`; it did not stop at the iteration cap. Its final classification changed 15,324 voxels from its prior state, with 0 lesion voxels on this no-lesion synthetic fixture.

Peak sampled RSS was **0.810 GB** under the 16 GB watchdog. A same-settings `--resume` verified all nine top-level stages and their output hashes; the resume invocation peaked at **0.243 GB** under a 4 GB watchdog. The audit validated 58 unique output paths by hash, 48 unique finite RAS NIfTI output paths across retained stage attempts, and three VTK surfaces tagged RAS. The top-level manifest names 17 final output records.

The 16 final scientific files match the earlier completed run on this same input/settings signature byte-for-byte. The Dura checkpoint manifest hash differs because it embeds paths under each run directory; it is metadata, not a changed scientific product.

## Evidence scope

This verifies fresh synthetic Python integration and restart behavior after the TGDM topology correction. It does **not** compare the complete MaCRUISEplus workflow to one Java/MATLAB execution and does not establish real-MRI validity. Selected component Java references remain separate evidence. Implementation code was unchanged during this run; implementation commit was `9eb84bc`. No full repository suite was run for this integration cycle; the latest documented suite predates that commit.

Detailed stage, convergence, RSS, output hash, and same-fixture comparison data are in the [machine-readable report](plus-after-tgdm-correction-2026-10-07.json).
