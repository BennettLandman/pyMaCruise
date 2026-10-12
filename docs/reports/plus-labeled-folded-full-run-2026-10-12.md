# Labeled folded synthetic PLUS full-run attempt

## Result

The current `plus-full-run` path completed input handling, orthogonal RAS
regridding, preprocessing, registration, registered recropping, first
membership initialization, and first thinning on the committed 65³ folded
synthetic fixture. It then stopped during the post-thinning GM-enabled
membership update because the source-faithful objective became nonfinite.
The manifest correctly records the run as failed. This is not a successful
end-to-end PLUS run and does not establish ribbon quality or source equivalence.

The fixture has labels 0/5/44/100 and deterministic tissue-specific T1/FLAIR
values, but it contains no lesion example. The bundled atlas descriptor
nevertheless enables lesion handling based on its topology declarations. The
first-iteration summary records two zero lesion-standard-deviation channels
and no zero standard-deviation tissue channels. VandyCRUISE's selected
`computeVariances` method assigns zero when a lesion denominator is empty;
Python's source-faithful IEEE handling then exposes the nonfinite objective at
the next membership pass. We will use a lesion-bearing synthetic case or
resolve the no-lesion source path from the Java wrapper before claiming this
checkpoint passes. The finite-zero-variance switch remains an explicit
comparison option and was not used here.

## Resources and artifacts

Peak sampled RSS was 3.378 GB under the 16 GB watchdog, so the failure was
algorithmic rather than a memory limit. Completed stages and the failure
manifest are retained under the ignored `build/full-plus-labeled-folded-current/`
directory. In particular, see:

- `plus-full-workflow.json`
- `lesion_workflow/first-iteration/first_iteration_summary.json`
- `lesion_workflow/registered-branch/lesion-registered-branch.json`

The input fixture remains at
`tests/data/synthetic/plus_labeled_folded_cortex/`. This run is synthetic
integration evidence only. It is distinct from the smaller CLI integration
test that completed but failed its independent ribbon-quality audit, and from
all Java source-comparison and real-human evidence.

## Reproduction

From the repository root:

```bash
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONNOUSERSITE=1 PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 1 --report-interval-seconds 60 -- \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 tests/data/synthetic/plus_labeled_folded_cortex/t1.nii.gz \
  --flair tests/data/synthetic/plus_labeled_folded_cortex/flair.nii.gz \
  --brain-mask tests/data/synthetic/plus_labeled_folded_cortex/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_labeled_folded_cortex/maseg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/full-plus-labeled-folded-current
```

## Next step

Trace how the V3 `run_maCRUISEplus` wrapper selects lesion-enabled behavior and
initial lesion memberships for a subject without a lesion, then create a
synthetic case that exercises that source path. Keep source-faithful NaN
behavior enabled until the Java path demonstrates otherwise.
