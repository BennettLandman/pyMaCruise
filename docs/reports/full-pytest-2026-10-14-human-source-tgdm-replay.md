# Full pytest: human source TGDM replay

- **Tested commit:** `1fa965f0520f5a8e42d79cdf75ca46bb0dec4bcf`
- **Result:** 1,577 passed, 4 skipped, 0 failed; 18 expected RAS/orientation warnings
- **Elapsed:** 849.60 seconds
- **Peak sampled RSS:** 10.549 GB, below the 16 GB watchdog; not terminated
- **Command:**

```bash
MPLCONFIGDIR=build/test-mpl-human-source-tgdm \
PYTHONNOUSERSITE=1 PYTHONPATH=src \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 1 --report-interval-seconds 30 -- \
  ./.venv/bin/python -m pytest -q
```

An initial collection-only launch used the standalone `pytest` executable, which omitted the repository root from the import path for `tests/test_prepare_plus_tgdm_java_inputs.py`. It stopped before running the suite. The successful full run above used `python -m pytest` and completed every test.
