# Full repository pytest — LesionTOADS repeated synthetic iteration

- Tested implementation commit: `c9c8b3abd00dd0759b004454afac553702bdea3d`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,167 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 877.24s (14:37)**.
- Peak sampled RSS: **14.369 GB**; the 24 GB watchdog did not terminate pytest.

This is repository regression evidence for P39bq's synthetic connected
iteration test, including registered atlas-shape sampling. It is not a
composite Java comparison, complete LesionTOADS equivalence, or real-image
validation. See the [P39bq report](lesion-registered-repeated-iteration-2026-10-04.md).
