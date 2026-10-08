# Full pytest validation: JIST connectivity-18 surface port

- Tested tree: `1bb0ab80fd9c26d2d139e32599ac27d42e472194` (plan/report commit); implementation: `194fb02c2e2b7b7bb9204f27a7487cadfb1ab380`.
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 0.2 --report-interval-seconds 60 -- ./.venv/bin/pytest -q`
- Result: **1,500 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 1,815.06 seconds**.
- Peak sampled RSS: **13.010 GB**; watchdog limit 16 GB; process completed normally.
- Output captured in the ignored build log `build/full-pytest-jist18-2026-10-08.txt`.

The full run includes the 256 per-case connectivity-table tests. It validates Python regressions at the tested tree; selected-JIST source equivalence is separately scoped in [the case-template comparison report](jist-connectivity18-case-template-port-2026-10-08.md), and whole-workflow MATLAB/JIST equivalence and real-MRI validation remain open.
