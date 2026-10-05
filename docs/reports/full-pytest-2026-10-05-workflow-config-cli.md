# Full PyMaCRUISE regression suite — 2026-10-05 config-driven workflow CLI

Implementation and package-verification commit: `308ffd52fa03428bf31dc0cc2dd37e63f0ea1951`.

Command:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,196 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
1,211.97s**. Peak sampled RSS: **13.827 GB**. The watchdog did not terminate
the run.

This is full Python regression-suite evidence for implementation commit
`308ffd5`. It includes synthetic base/plus workflows, including complete
plus run launched through the JSON-config CLI and resumed, and the focused
config-dispatch behavior described in
[the P42a report](workflow-config-cli-2026-10-05.md). It is not source-wide
Java equivalence or real-MRI validation.
