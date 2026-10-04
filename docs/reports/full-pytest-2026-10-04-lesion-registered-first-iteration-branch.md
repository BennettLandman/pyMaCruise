# Full repository pytest — LesionTOADS registered first iteration branch

- Tested implementation commit: `6488c76d96d7342530a38b23f39f389b471622f8`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,165 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 876.76s (14:36)**.
- Peak sampled RSS: **13.695 GB**; the 24 GB watchdog did not terminate pytest.

This is repository regression evidence for P39bo's synthetic composition. It
is not full LesionTOADS equivalence or real-image validation. See the [P39bo
report](lesion-registered-first-iteration-branch-2026-10-04.md).
