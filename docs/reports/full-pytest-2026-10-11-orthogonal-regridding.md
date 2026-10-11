# Full-suite validation: optional orthogonal RAS regridding

Tested implementation tree: `cabb52aa` (includes implementation commit
`689465d` and its focused-validation report commit).

## Result

```text
1554 passed, 4 skipped, 0 failed, 11 expected RAS-orientation warnings
858.73 seconds
```

The suite ran under `scripts/run_with_rss_watchdog.py` with a **16 GB** ceiling.
Peak sampled RSS was **9.443 GB**; the watchdog did not terminate the test
process. The expected warnings are from tests that deliberately provide LAS or
LPS images and verify RAS normalization.

## Reproduction

From `pymacruise/`:

```sh
MPLCONFIGDIR=build/full-suite-mpl-regrid PYTHONNOUSERSITE=1 PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 120 -- \
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q
```

The suite includes the oblique synthetic base full-workflow/restart case, plus
regridding resource/hash and opt-out checks, shared image-I/O tests, and the
existing source-comparison and workflow regression suites. This is automated
software/synthetic validation. It does not establish source equivalence or
real-human output acceptance.
