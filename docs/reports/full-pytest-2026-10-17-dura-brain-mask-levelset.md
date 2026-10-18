# Full pytest after BASE Dura brain-mask level-set routing — 2026-10-17

## Result

The full test suite passed on implementation commit `8438185d49fcb8ebd49626dc249d763b89a02ffc` (`Preserve Dura brain mask level set through base TGDM`).

- **1,683 passed, 4 skipped, 0 failed**
- **30 warnings**, emitted by the expected RAS reorientation paths
- **Runtime:** 1,000.52 seconds (16 minutes 40.52 seconds)
- **Peak sampled RSS:** 9.435 GB under a 16 GB watchdog; it did not terminate the run
- **Focused integration selection:** 43 passed in 44.48 seconds
- **Command:**

```bash
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- \
  ./.venv/bin/pytest -q --tb=short
```

The captured runner log is `build/full-pytest-2026-10-17-dura-brain-mask-levelset.log` (ignored build output). This is Python regression evidence on the stated commit. It does not establish a fresh full original MaCRUISE reference run.
