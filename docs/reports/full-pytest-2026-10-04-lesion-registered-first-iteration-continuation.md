# Full repository pytest — LesionTOADS continuing first iteration

- Tested implementation commit: `629ec61127b541153c49bb766b16a2b40f9d0466`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,166 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 850.37s (14:10)**.
- Peak sampled RSS: **14.717 GB**; the 24 GB watchdog did not terminate pytest.

This is repository regression evidence for the synthetic continuing-branch
composition. It does not establish a converged multi-iteration run, complete
LesionTOADS equivalence, or real-image validation. See the [P39bp report](lesion-registered-first-iteration-continuation-2026-10-04.md).
