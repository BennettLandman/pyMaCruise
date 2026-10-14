# Full pytest after MaCRUISE V3 Dura parameter correction — 2026-10-14

The complete suite passed against the implementation and test tree at commit
`9ee3ee92965a7b71f8107086dc51334f9ffd6935`.

```text
1,577 passed, 4 skipped, 0 failed, 18 expected orientation warnings
719.49 seconds (11 minutes 59 seconds)
Peak sampled RSS: 15.100 GB under a 16 GB watchdog
```

Command, from the repository root:

```sh
PYTHONPATH=.:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/pytest -q
```

The suite includes the corrected synthetic Dura wrapper fixtures, both WM
0.05/GM 0.02 topology parameters in CLI provenance, and restart/resume command
serialization. The suite passed with roughly 1.1 GB of sampled headroom under
the configured watchdog. This is code/test validation; the refreshed human
BASE workflow and its second TOADS pass are still pending.
