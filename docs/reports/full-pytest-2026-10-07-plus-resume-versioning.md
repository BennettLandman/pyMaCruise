# Full pytest after plus checkpoint versioning — 2026-10-07

## Command and result

Tested commit: `c17d2c4` (includes implementation `caf876d`, which adds implementation-versioned plus stage generations, and its regression test).

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/pytest -q --tb=short
```

Result: **1,230 passed, 0 failed, 0 skipped**, with **11 expected RAS-orientation warnings**, in **1,607.32s (26m 47s)**. The watchdog sampled a peak RSS of **13.364 GB**, below its 16 GB limit. The process exited normally; the watchdog did not terminate it.

## Scope

This full-suite run covers the plus checkpoint-generation change and the existing integration regressions. The focused generation-change test passed 1/1 in 102.89s at 0.582 GB; the existing restart/corruption-repair test passed 1/1 in 79.55s at 0.579 GB on implementation commit `caf876d`.

No new VandyCRUISE or other algorithm source methods were mapped in this implementation cycle. The method-based selected TOADS estimate remains **3,511–4,911 of 13,651 source lines (about 26–36%)**; this is not runtime Python coverage and is unchanged by the orchestration/cache correction. The run demonstrates Python regression-suite success, not source equivalence or real-MRI validation.
