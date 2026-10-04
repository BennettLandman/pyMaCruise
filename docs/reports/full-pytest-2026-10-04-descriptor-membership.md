# Full pytest after LesionTOADS descriptor membership support

Command:

```text
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

Tested commit: `11cb482c17282589ebe8777406644670da69f189`.

Result: **1,149 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
785.28 seconds (13:05)**. The watchdog exited 0, sampled peak RSS of
**15.498 GB**, and did not terminate the child process (24 GB cap).

The suite ran synthetic test data and did not execute the supplied human MRI.
This validates regressions after relaxing the shared membership-kernel
validation for the bundled 13-class descriptor, together with its new
bitwise one-voxel Java formula test. It does not validate a connected
LesionTOADS iteration or the full workflow. No new selected base-TOADS source
lines were counted; its coarse estimate remains 3,464–4,864 of 13,651
(25–36%).
