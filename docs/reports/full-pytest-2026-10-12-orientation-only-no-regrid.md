# Full-suite validation: orientation-only no-regrid behavior

Implementation commit: `dcf1a6cf` (`Clarify orientation-only no-regrid mode`).

## Command and result

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- \
  .venv/bin/pytest -q
```

Result: **1,563 passed, 4 skipped, 0 failed** in **852.83 seconds**. Pytest
reported 18 expected input-orientation warnings, including the new warnings
that identify oblique direction-cosine removal in the no-regrid synthetic
cases. The 16 GB watchdog sampled a peak RSS of **15.493 GB** and did not
terminate the run.

## Scope

This is the first full-suite run covering the clarified opt-out semantics:
closest-RAS axis permutations/flips, no interpolation or voxel-value change,
orthogonal RAS output headers, center-anchored origin, and explicit warnings
when oblique directions are discarded. It also rechecks default-on resampling,
the BASE/PLUS workflow tests, resume behavior, and all other package tests.
The full suite validates the synthetic/software contracts; the replacement
human BASE run using the clarified option is tracked separately in P46-H2.
