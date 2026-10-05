# Full pytest — LesionTOADS Dura output checkpoint

- Tested implementation commit: `2b1c92ceddc7400016f1d0ccb70a32b2d9740419`
- Command: `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,193 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Duration: **1,000.66 seconds**
- Peak sampled RSS: **14.374 GB**; configured limit **24 GB**; terminated for RSS limit: **false**
- Classification: repository regression evidence. The P41r focused synthetic test and row-specific report describe its narrower workflow evidence; the full suite does not establish whole Java equivalence or real-image validation.
