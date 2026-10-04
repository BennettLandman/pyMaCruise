# Full repository pytest — LesionTOADS repeated synthetic iteration

- Tested implementation commit: `d962c50be6b0ecdef0f8d79415117854edced255`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,167 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 882.79s (14:42)**.
- Peak sampled RSS: **14.513 GB**; the 24 GB watchdog did not terminate pytest.

This is repository regression evidence for P39bq's synthetic connected
iteration test. It is not a composite Java comparison, complete LesionTOADS
equivalence, or real-image validation. See the [P39bq report](lesion-registered-repeated-iteration-2026-10-04.md).
