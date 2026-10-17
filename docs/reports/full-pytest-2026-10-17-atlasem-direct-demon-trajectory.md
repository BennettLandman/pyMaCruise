# Full test suite — AtlasEM optimizer trajectory regression

- **Tested implementation commit:** `9245b50689ff17e317aa22ff2058ed65edc01d09`
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-17-atlasem-demon-trajectory .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,651 passed, 4 skipped, 0 failed; 21 expected orientation warnings.
- **Elapsed:** 774.40 s (12 min 54.40 s).
- **Peak sampled RSS:** 10.037 GB; the 16 GB watchdog did not trigger.

The suite includes direct execution of selected Java
`DemonToadDeformableAtlas` and exact Python transform comparison for synthetic
single-scale and pyramid cases. The pyramid test includes nonzero updates at
iteration caps 1, 2, and 5. Selected atlas-resource loading and complete
AtlasEM/PLUS behavior remain outside this evidence.
