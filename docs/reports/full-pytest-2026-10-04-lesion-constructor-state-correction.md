# Full pytest evaluation — LesionTOADS constructor-state correction

The repository-wide suite passed against implementation commit `44c2916`:

```text
1,159 passed, 11 warnings in 862.74s (0:14:22)
```

The warnings are expected RAS-reorientation notices from orientation-contract
tests. The suite used synthetic fixtures; it did not run the supplied human
MRI. The 24 GB RSS watchdog completed normally with peak sampled RSS of
14.140 GB.

Command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q
```

This suite includes the corrected LesionTOADS constructor-state reference and
its first-thinning integration test. It does not establish complete
LesionTOADS or MaCRUISEplus equivalence.

The selected TOADS source-line estimate remains **3,464–4,864 of 13,651
(about 25–36%)**, with zero newly counted selected base-TOADS method lines in
this cycle. This is coarse method-level accounting, not runtime Python code
coverage or workflow completion.
