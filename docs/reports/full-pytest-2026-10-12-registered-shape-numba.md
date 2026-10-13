# Full pytest evaluation — registered-shape Numba backend

- **Implementation commit:** `645b56b2` (`Accelerate registered atlas prior sampling`)
- **Documentation/status commit:** `8c5fcd35` (`Record registered-loop quality plateau`)
- **Tested HEAD:** `8c5fcd35b9f872c8d5e41a5683479d5ff19d944f`
- **Result:** **1,568 passed, 4 skipped, 0 failed**; 18 expected orientation warnings.
- **Elapsed:** 1,592.23 seconds (26:32).
- **Peak sampled RSS:** 12.990 GB under a 16 GB watchdog; the watchdog did not terminate the run.
- **Environment:** Python virtual environment; user site disabled; OMP, OpenBLAS, and MKL thread counts set to 1; 16 GB watchdog.

## Reproduction

```sh
env PYTHONPATH=src PYTHONNOUSERSITE=1 MPLCONFIGDIR="$PWD/build/matplotlib-config" \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 2 --report-interval-seconds 30 -- \
  .venv/bin/pytest -q --tb=short
```

The run validates the complete existing test suite, including the new Numba-versus-scalar exact-output tests and existing registered-shape Java harness comparison. This is test-suite evidence; it does not turn the atlas-derived synthetic PLUS run into independent ground truth or resolve the outstanding ribbon-quality failures. The prior suite result of 1,563 passed on `dcf1a6cf` is superseded for the current code state.
