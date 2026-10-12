# Full suite after source-faithful LesionTOADS NaN continuation

## Result

On implementation/documentation commit `30c46def`, the complete repository
suite passed **1,559 tests, 4 skipped, 0 failed**, with 11 expected RAS
orientation warnings in **824.78 seconds**. Peak sampled RSS was **11.880 GB**
under a 16 GB watchdog; the watchdog did not terminate the run.

Command from the repository root:

```bash
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONNOUSERSITE=1 PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 1 --report-interval-seconds 60 -- \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q
```

Focused source-faithful LesionTOADS tests passed **19 tests in 39.61 seconds**
on implementation commit `9d7c7c31`. In addition, the public `plus-full-run`
completed all nine stages on the labeled 65³ synthetic case using default-on
RAS regridding and a one-pass cap. That separate run peaked at **5.433 GB** and
passed all 51 manifest hash checks and all 46 RAS NIfTI orientation checks;
its ribbon-quality audit failed. The combined results remain synthetic
integration evidence, not human PLUS acceptance or whole-workflow Java
equivalence. See [the public run report](plus-labeled-folded-full-run-2026-10-12.md).

## Source-line estimate

This cycle maps **0 additional selected TOADS source-method lines** to the
coarse method estimate. It remains **3,511–4,911 of 13,651 lines (26–36%)**.
This is source-method mapping, not Python runtime line coverage or workflow
completion.
