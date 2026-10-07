# Full pytest — plus TGDM uncrop exterior correction

- Tested commit: `7ae1dee` (`Compare plus meshes with selected JIST source`)
- Code correction: `2f863b4` (`Fix plus TGDM uncrop exterior`)
- Command:

  ```bash
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    PYTHONNOUSERSITE=1 PYTHONPATH=src \
    MPLCONFIGDIR=/private/tmp/pymacruise-full-suite-mpl \
    .venv/bin/python scripts/run_with_rss_watchdog.py \
      --limit-gb 24 --interval-seconds 1 --report-interval-seconds 60 -- \
      .venv/bin/pytest -q --tb=short
  ```

- Result: **1,220 passed, 0 failed, 0 skipped, 11 expected RAS warnings**
- Duration: **1,463.83 s (24m 23.83s)**
- Peak sampled RSS: **12.448 GB** under the 24 GB watchdog; child exit 0,
  watchdog did not terminate the run.

This suite includes the 20-valued MaTGDMcrop uncrop exterior correction and
all repository tests on the tested commit. It validates Python regressions;
it does not establish whole-workflow MATLAB/JIST equivalence or real-MRI
validation.
