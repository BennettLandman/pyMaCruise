# Full repository pytest — LesionTOADS third synthetic iteration

- Tested implementation commit: `6def02db2281f804efbcc40d506bf3fd7ebcb7c5`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,167 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 903.35s (15:03)**.
- Peak sampled RSS: **14.892 GB**; the 24 GB watchdog did not terminate pytest.

This repository regression run includes the three-update synthetic
LesionTOADS loop test. It is not source-equivalence evidence for the composed
workflow, full LesionTOADS validation, or real-image validation. See the
[third-iteration report](lesion-registered-third-iteration-2026-10-04.md).
