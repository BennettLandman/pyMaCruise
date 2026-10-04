# Full pytest after P41a plus integration

Command:

```text
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Tested implementation commit: `677975369a8a11e7776717eed51a44dcee517c69`.

Result: **1,154 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
834.21 seconds (13:54)**. Peak sampled RSS was **15.381 GB** under the 24 GB
watchdog; it exited successfully and did not terminate the test process.

This regression run includes P41a, which connects synthetic LesionTOADS Dura
channels through transfer/crop, selected topology and pre-TGDM stages, nested
TGDM, uncropping, surface extraction, and thickness. Atlas and TOADS images are
committed synthetic inputs, not outputs from a complete LesionTOADS execution.
The run does not validate the full plus workflow, FLAIR/atlas registration,
restart/manifest behavior, Java whole-workflow equivalence, or real data.
No selected base-TOADS source lines were newly counted; its coarse estimate
remains 3,464–4,864 of 13,651 (25–36%).
