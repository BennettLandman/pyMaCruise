# Full repository pytest — LesionTOADS registered iteration

- Tested implementation commit: `73056475da515127c43aa64ade613a15cddd3197`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,159 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 862.69s (14:22)**.
- Peak sampled RSS: **14.584 GB**; the 24 GB watchdog did not terminate pytest.

This is repository regression evidence for the implementation commit, not
whole-workflow source equivalence or real-image validation. The focused
13-class registered iteration comparison is in the
[slice report](lesion-descriptor-registered-iteration-2026-10-04.md).
