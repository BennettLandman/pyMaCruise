# Full test suite after N4 ANTs CLI default correction

## Result

On implementation commit `f77e5b54` (`Match N4 intensity rescaling to ANTs CLI`), the full repository suite passed: **1,548 passed, 4 skipped, 0 failed**, with **11 expected RAS-orientation warnings** in **830.03 seconds**. Peak sampled RSS was **10.487 GB** under a **16 GB** watchdog; the watchdog did not terminate the run.

This validates the repository tests after the shared `N4Config` default changed to intensity rescaling enabled. It is regression evidence, not source equivalence or a human-workflow acceptance pass. The downloaded source N4 output still differs from Python by 11.75% brain-region relative MAE with rescaling enabled; see the [same-scan report](in-vivo-source-output-comparison-2026-10-11.md).

## Reproduction

From `pymacruise/`:

```sh
MPLCONFIGDIR=build/full-suite-mpl PYTHONNOUSERSITE=1 PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 120 -- \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q
```
