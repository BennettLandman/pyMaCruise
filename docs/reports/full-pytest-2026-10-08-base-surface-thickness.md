# Full pytest after base surface/thickness integration — 2026-10-08

## Command and result

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/python -m pytest -q
```

Tested commit: `5a412b7`.

Result: **1,243 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 1,743.86s** (29m 3.86s). Peak sampled RSS was **14.961 GB**. The 16 GB watchdog did not terminate the run.

This suite includes the seven-stage small-fixture base workflow/restart test, including finite RAS thickness volume, RAS-mm surfaces with closest-point fields, checkpoint resume, and downstream invalidation. The full-size surface/thickness component output is separately documented in [the brain-resolution report](base-surface-thickness-brain-resolution-2026-10-08.md). It does not establish a full-size seven-stage manifest, whole-workflow Java/MATLAB equivalence, or real-MRI validation. No selected Java algorithm source lines were added in this cycle; the coarse selected TOADS source-method estimate remains 3,511–4,911 of 13,651 lines (about 26–36%).
