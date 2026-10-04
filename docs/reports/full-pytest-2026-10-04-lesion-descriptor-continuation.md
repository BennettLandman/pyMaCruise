# Full repository pytest — LesionTOADS descriptor continuation

- Tested implementation commit: `5f7fc042ba8394435648d9d3b6f6e701d1d5f4df`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,159 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 876.20s (14:36)**.
- Peak sampled RSS: **14.697 GB**; the 24 GB watchdog did not terminate pytest.

This is repository regression evidence on the implementation commit. It does
not establish full MaCRUISE or LesionTOADS source equivalence, nor real-image
validation. The focused 13-class continuation comparison is documented in
[the slice report](lesion-descriptor-continuing-iteration-2026-10-04.md).
