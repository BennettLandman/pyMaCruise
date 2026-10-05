# Full pytest — plus full workflow integration

- Tested implementation commit: `7dc8a67b10489bbf2ba2e7d743c106366c208ae4`
- Command: `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,193 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Duration: **995.58 seconds**
- Peak sampled RSS: **14.609 GB**; configured limit **24 GB**; terminated for RSS limit: **false**
- Classification: repository regression evidence. It does not establish full MaCRUISE/plus Java equivalence or real-image validation.
