# Full test suite — AtlasEM multiple registered classes

- **Tested implementation commit:** `202939fd224e87298cd1f6cd931292f3e454f38c`
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-17-atlasem-multiple-registered .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,653 passed, 4 skipped, 0 failed; 21 expected orientation warnings.
- **Elapsed:** 769.95 s (12 min 49.95 s).
- **Peak sampled RSS:** 12.158 GB; the 16 GB watchdog did not trigger.

This suite includes direct selected-Java and Python pyramid optimizer
comparisons for three-class inputs with one or two non-background registered
classes. The tests verify bit-identical, nonzero final rigid transforms.
Selected descriptor resources and complete AtlasEM/PLUS behavior remain
outside this evidence.
