# Full repository pytest — LesionTOADS repeated synthetic iteration

- Tested implementation commit: `37b857966fda75c04772541d0c8170b1db1cde68`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,167 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 868.23s (14:28)**.
- Peak sampled RSS: **13.738 GB**; the 24 GB watchdog did not terminate pytest.

This is repository regression evidence for P39bq's synthetic connected
iteration test, including refreshed shape maps, second registration, and the
terminal classification branch. It is not a composite Java comparison,
complete LesionTOADS equivalence, or real-image validation. See the [P39bq
report](lesion-registered-repeated-iteration-2026-10-04.md).
