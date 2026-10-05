# Full PyMaCRUISE regression suite — 2026-10-05 adaptive LesionTOADS fields

Tested implementation commit: `2b90545a52ace6a3d518925a380c9c34ccac88a0`.

Command:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,200 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
1,239.07s**. Peak sampled RSS: **12.993 GB**. The watchdog did not terminate
the run.

This suite includes the fixed-field standard-13 membership path and the
opt-in adaptive correction-field refit in the registered LesionTOADS loop.
The adaptive feature is an in-memory API capability; this run does not add
checkpoint persistence or CLI/config field options and is not a composite
Java/source-equivalence result. No human MRI was used.
