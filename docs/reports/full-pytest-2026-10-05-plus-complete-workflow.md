# Full pytest — raw-input MaCRUISEplus workflow

- Tested implementation commit: `3051b821ed93163f750b4e2b81fb98b0bad22a52`
- Command: `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,193 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Duration: **1,079.02 seconds**
- Peak sampled RSS: **15.100 GB**; configured limit **24 GB**; terminated for RSS limit: **false**
- Classification: repository regression evidence. It does not establish whole MaCRUISEplus Java equivalence or real-image validation.
