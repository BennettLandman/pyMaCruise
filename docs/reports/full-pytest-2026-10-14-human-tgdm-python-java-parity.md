# Full pytest: human TGDM Python/Java same-input comparison

The full repository suite passed on tested tree commit `8fd85e9467fa3d8be7df6b1a2e21a0a26e5f28c5`:

- **1,577 passed, 4 skipped, 0 failed**
- **18 expected orientation warnings**
- **762.82 seconds** (12 minutes 42.82 seconds)
- **13.731 GB peak sampled RSS**, under a 16 GB watchdog; not terminated

Run from the `pymacruise` repository root:

```bash
MPLCONFIGDIR=build/test-mpl-human-source-tgdm-python-parity \
PYTHONNOUSERSITE=1 PYTHONPATH=src \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 30 -- \
  ./.venv/bin/python -m pytest -q
```

The 18 warnings are expected RAS reorientation/orientation-only notices from their dedicated tests. This suite validates the repository after adding the direct same-input human TGDM replay; the replay's numerical source comparison is reported separately in [the TGDM report](in-vivo-human-source-tgdm-replay-2026-10-14.md).
