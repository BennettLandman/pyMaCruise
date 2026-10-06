# Natural-convergence synthetic MaCRUISEplus workflow — 2026-10-06

## Result

The default raw-input `plus-full-run` path now iterates the registered
LesionTOADS loop until its source tolerance is reached or the configured
iteration ceiling is reached. The default P41k/P41m tolerance is `1e-4`; the
complete-run pass ceiling is 99. The explicitly named
`plus-lesion-registered-loop` remains a bounded checkpoint command with a
one-pass default.

The committed 36×45×34 asymmetric folded RAS fixture ran all nine workflow
stages. LesionTOADS stopped naturally at iteration 5 with energy ratio
`8.241410796472337e-06` (below `1e-4`), and wrote terminal classification and
lesion-subtype images. The test verifies stage/output hashes, finite RAS
NIfTIs, native T1-grid products, final surface/thickness outputs, and resume
with the exact original settings. Focused run peak sampled RSS was 8.798 GB
under the 16 GB watchdog. This is synthetic Python integration evidence, not
Java/MATLAB whole-workflow equivalence, numerical acceptance against the
source, brain-resolution plus validation, or real-MRI validation.

## Focused validation

Command:

```bash
env MPLCONFIGDIR=/private/tmp PYTHONPATH=src \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pytest -q --tb=short \
  tests/test_plus_full_workflow.py::test_plus_full_run_completes_on_larger_asymmetric_folded_ras_grid
```

Result on implementation commit `bb40acc0d0b618f62d035bc3eb32f8a313567e50`:
**1 passed, 0 failed, 0 skipped in 469.67 seconds**, peak sampled RSS
**8.798 GB**. The test reaches convergence at iteration 5 and verifies
restartability on the same complete-run defaults.

## Full repository regression

Command:

```bash
env MPLCONFIGDIR=/private/tmp PYTHONPATH=src \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  ./.venv/bin/pytest -q --tb=short
```

Full suite on tested implementation/documentation commit
`09a7bd034f1d32c336525a9d50002af897f79140`: **1,219 passed, 0 failed,
0 skipped, 11 expected RAS warnings in 1,699.06 seconds**. Peak sampled RSS was
**13.052 GB** under the 24 GB watchdog.

This workflow integration adds zero selected VandyCRUISE TOADS source lines.
After the full-suite evaluation the method-level selected TOADS source
representation estimate remains **3,511–4,911 of 13,651 lines (about
26–36%)**; this coarse source estimate is not Python test coverage or
workflow completion.
