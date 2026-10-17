# Full test suite — AtlasEM selected-class registration

- **Tested implementation commit:** `f727c51431e19706433f262b8838f1be74d58fe4`
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-17-atlasem-class-selection .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,652 passed, 4 skipped, 0 failed; 21 expected orientation warnings.
- **Elapsed:** 743.94 s (12 min 23.94 s).
- **Peak sampled RSS:** 13.114 GB; the 16 GB watchdog did not trigger.

The suite includes direct execution of selected Java
`DemonToadDeformableAtlas` and exact Python transform comparison on a
three-class pyramid case with one selected non-background class. The
selected `ibsr8obj` resources and complete AtlasEM/PLUS behavior remain
outside this evidence.
