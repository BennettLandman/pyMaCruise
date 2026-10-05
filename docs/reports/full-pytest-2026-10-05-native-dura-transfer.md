# Full pytest — native-space Dura transfer

- Tested implementation commit: `2c32fd89ae62d75ffa3e445e8f7079e59b620d10`
- Command: `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,193 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Duration: **1,026.32 seconds**
- Peak sampled RSS: **14.479 GB**; configured limit **24 GB**; terminated for RSS limit: **false**
- Classification: repository regression evidence. It does not establish JIST arbitrary-grid transfer equivalence, full MaCRUISEplus source equivalence, or real-image validation.
