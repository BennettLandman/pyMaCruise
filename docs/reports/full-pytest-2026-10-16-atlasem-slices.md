# Full test suite — AtlasEM fixed-input source slices

- **Tested tree commit:** `2aa83740` (`Port AtlasEM Dura outputs and MIPAV hole fill`)
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-16-atlasem-slices .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,598 passed, 4 skipped, 0 failed; 18 expected orientation warnings.
- **Elapsed:** 780.86 s (13 min 00.86 s).
- **Peak sampled RSS:** 11.981 GB; the 16 GB watchdog did not trigger.
- **Selected TOADS source-method estimate:** unchanged at 3,511–4,911 of 13,651 Java lines (26–36%). This checkpoint added no TOADS methods; it added fixed-input AtlasEM slices and MIPAV 2.5D hole-fill behavior.

The seven additional focused test cases cover AtlasEM fixed-input statistics,
one single-channel Atlas-prior membership update, AtlasEM Dura-input products,
and MIPAV fill-holes boundary/flood semantics. They use committed synthetic
arrays and pinned Java source excerpts. This run does not execute or validate
the missing AtlasEM atlas assets, deformable registration, complete EM stage,
or pre-LesionTOADS Dura evolution. See [P44b source-stage audit](plus-source-stage-coverage-audit-2026-10-16.md).
