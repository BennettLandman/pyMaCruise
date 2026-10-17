# Full test suite — AtlasEM direct Demon single-scale and pyramid optimizers

- **Tested code/test commit:** `f179290210fb2f0c9cdb193833631dd473081799`
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-16-atlasem-demon-pyramid .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,650 passed, 4 skipped, 0 failed; 21 expected orientation warnings.
- **Elapsed:** 719.57 s (11 min 59.57 s).
- **Peak sampled RSS:** 14.822 GB; the 16 GB watchdog did not trigger.

The suite includes direct selected-Java comparisons for the single-scale
post-recrop and initial pyramid AtlasEM optimizers. The six final float32
transform parameters match bit for bit for the tested synthetic cases.
This result does not supply the unavailable selected `ibsr8obj`
descriptor/payload or validate complete AtlasEM/PLUS behavior.
