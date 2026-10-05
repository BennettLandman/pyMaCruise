# Full PyMaCRUISE regression suite — 2026-10-05 LesionTOADS field checkpoints

Tested implementation commit: `cd6d3f2cd76b0583bad67466a22c34085c8ffe24`.

Command:

```bash
PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Result: **1,201 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
1,205.29s**. Peak sampled RSS: **10.986 GB**. The watchdog did not terminate
the run.

The implementation carries optional channel-last correction fields and
adaptive-fit degree through registered LesionTOADS checkpoint files, resumes
them in later loop segments, and exposes the controls on both the checkpoint
CLI and restartable LesionTOADS loop workflow. Tests use synthetic NIfTI data.
This is not whole-pipeline Java equivalence or human-MRI validation.
