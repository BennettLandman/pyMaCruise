# Full PyMaCRUISE suite — oblique no-regrid support

## Result

On implementation commit `452cc573`, the full repository suite passed
**1,561 tests, 4 skipped, 0 failed**, with 11 expected RAS orientation
warnings in **829.51 seconds**. Peak sampled RSS was **12.018 GB** under a
16 GB watchdog; the watchdog did not terminate the run.

Command from the repository root:

```bash
PYTHONPATH=src PYTHONNOUSERSITE=1 \
  MPLCONFIGDIR=build/test-mpl-full-oblique-optout \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 45 -- \
  .venv/bin/pytest -q
```

The focused geometry/workflow selection passed **19 tests in 0.94 seconds**
with peak sampled RSS **0.125 GB**. It verifies default-on resampling,
oblique-grid preservation under `--no-regrid-inputs`, manifest/resource
provenance, TOADS signed voxel spacings, and shear rejection. See the updated
[regridding report](orthogonal-ras-regridding-focused-2026-10-11.md).

## Source-method coverage

This geometry-path change maps **0 additional selected TOADS source-method
lines**. The coarse estimate remains **3,511–4,911 of 13,651 lines (26–36%)**;
it is method mapping, not runtime Python coverage or workflow completion.

## Evidence limits

The suite validates the current Python regressions and synthetic no-regrid
workflow boundaries. It does not compare an oblique no-regrid human workflow
against a complete historical Java/MATLAB run, resolve the existing same-scan
human mismatch, or validate the PLUS workflow on human MRI.
