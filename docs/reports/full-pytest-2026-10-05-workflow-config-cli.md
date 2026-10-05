# Full PyMaCRUISE regression suite — 2026-10-05 config-driven workflow CLI

Implementation commit: `4b86a40f3247e09c3a19949f7e0ddd2485dc0247`.

Command:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,195 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
1,083.83s**. Peak sampled RSS: **15.102 GB**. The watchdog did not terminate
the run.

This is full Python regression-suite evidence for implementation commit
`4b86a40`. It includes synthetic base/plus workflows and the focused
config-dispatch behavior described in
[the P42a report](workflow-config-cli-2026-10-05.md). It is not source-wide
Java equivalence or real-MRI validation.
