# Full pytest after mask-aware human PLUS TGDM correction — 2026-10-17

## Result

The full test suite passed on implementation commit `50ad2ea1b25ecde845bc8c399813c9da3297b419` (`Route brain mask through PLUS TGDM source comparison`).

- **1,681 passed, 4 skipped, 0 failed**
- **30 warnings**, all emitted by expected input reorientation paths that report non-RAS inputs
- **Runtime:** 991.89 seconds (16 minutes 31.89 seconds)
- **Peak sampled RSS:** 7.758 GB under a 16 GB watchdog; the watchdog did not terminate the run
- **Command:**

```bash
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- \
  ./.venv/bin/pytest -q --tb=short
```

The captured log is `build/full-pytest-2026-10-17-mask-aware-tgdm.log` (ignored build output). This suite validates regression coverage for the current Python checkout. The human Java comparison remains selected-stage, held-fixed TGDM evidence; it is not a complete independent source workflow comparison. The matching mask input is the Python adapter's unit signed bound, and its equivalence to the upstream source `Brain Mask Level Set` remains unverified.
