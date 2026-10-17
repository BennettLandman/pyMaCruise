# Full test suite — PLUS AtlasEM/Dura stage integration

- **Tested implementation commit:** `51e87d3b15b6792313edb1a3b1783880b276778f`
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-17-plus-atlasem-dura-integration .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,673 passed, 4 skipped, 0 failed; 21 expected orientation warnings.
- **Elapsed:** 665.27 s (11 min 5.27 s).
- **Peak sampled RSS:** 15.388 GB; the 16 GB watchdog did not trigger.

The suite includes generated AtlasEM/Dura handoff, resume, and parent-stage
ordering tests. It validates synthetic integration and legacy behavior, not
the unavailable selected `ibsr8obj` assets or whole PLUS source equivalence.
