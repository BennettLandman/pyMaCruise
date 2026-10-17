# Full regression evaluation — TOADS compiled execution — 2026-10-17

- Result: **1,670 passed, 4 skipped, 0 failed**, with 21 expected orientation warnings.
- Duration: **721.40 s** (12 minutes 1.40 seconds).
- Peak sampled child RSS: **8.475 GB** under a 16 GB watchdog; no RSS termination.
- Production numerical kernels match commit `6c3ad57a` (compiled implementations in `f2158a79`, source/reference regressions in `6c419308`, NaN behavior in `6c3ad57a`). Subsequent edits to these production wrappers add documentation only.

The suite includes the existing selected Java field, growing, relation, registration, and nineteen-iteration membership/pose/classification regressions, along with BASE/PLUS integration checks. This is regression evidence and does not close complete external BASE/PLUS source equivalence.

The full suite collected the initial 17 new compiled/reference tests. The additional zero-threshold NaN-priority regression was added after collection and passed in the separate focused run of **18 tests in 4.38 s**. Thus the extra edge case is covered without claiming it was part of the earlier collection.

The initial full-suite command used `PYTHONPATH=src` and stopped during collection because `scripts.prepare_plus_tgdm_java_inputs` could not be imported. The corrected command completed successfully:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 MPLCONFIGDIR=build/test-mpl \
  .venv/bin/python -u scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 30 -- \
  .venv/bin/pytest -q --tb=short
```

Raw output is ignored at `build/toads-performance-2026-10-17/full-suite-correct-path.log`. A fresh human TOADS replay ran concurrently under its own 12 GB watchdog, so the elapsed suite duration is not an isolated performance benchmark. TOADS performance is reported separately in [the optimization evaluation](toads-compiled-performance-2026-10-17.md). The source-method line estimate is unchanged because these kernels accelerate methods already mapped in the plan.
