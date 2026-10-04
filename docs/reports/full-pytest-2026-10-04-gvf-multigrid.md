# Full pytest evaluation after VandyCRUISE GVF multigrid port

## Result

The full suite passed on implementation commit
`ea31121f2366408666f28ecd11c2008f83c81859`:

```text
1141 passed, 11 warnings in 801.80s (0:13:21)
```

There were no failed or skipped tests. The 11 warnings are expected RAS
reorientation notices from orientation-contract tests. No real MRI was used.

Command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

The watchdog completed normally (`child_exit=0`, not terminated for memory)
and sampled peak RSS of **15.760 GB**. This is below the configured 16 GB
limit, but close enough that later test cycles should favor focused suites
unless a shared algorithm contract or major workflow milestone warrants a
full run. The peak is also below the machine's 40 GB maximum.

## Scope of evidence

This is repository-wide regression evidence for the Numba port of the selected
VandyCRUISE Poisson3d multigrid V-cycle. It does not add mapped selected TOADS
Java source lines, so the coarse source-line estimate remains
**3,464–4,864 of 13,651 (25–36%)**, with zero source-line delta this cycle.
The separate full-grid Java comparison and the restartable base-run evidence
remain in [the source-atlas report](toads-source-vandy-atlas-comparison-2026-10-03.md).
Synthetic suite success does not establish TOADS convergence, whole-workflow
Java equivalence, or real-data validation.
