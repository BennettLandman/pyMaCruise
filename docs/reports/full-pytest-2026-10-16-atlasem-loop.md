# Full test suite — AtlasEM fixed-prior iteration loop

- **Tested commit:** `cc82ba72` (`Port AtlasEM fixed-prior iteration loop`)
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-16-atlasem-loop .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,618 passed, 4 skipped, 0 failed; 18 expected orientation warnings.
- **Elapsed:** 777.11 s (12 min 57.11 s).
- **Peak sampled RSS:** 5.120 GB; the 16 GB watchdog did not trigger.
- **Selected TOADS source-method estimate:** unchanged at 3,511–4,911 of 13,651 Java lines (26–36%). AtlasEM work adds no methods to the selected TOADS source estimate. This is a source-method mapping estimate, not runtime coverage or workflow completion.

The tested code includes source-excerpt validated synthetic AtlasEM fixed-input
statistics, preprocessing, membership updates, fixed-prior loop behavior, and
output mappings. It does not validate the absent `ibsr8obj` payload, atlas
registration, the complete AtlasEM stage, pre-LesionTOADS Dura evolution, or
whole PLUS equivalence. See [P44b source-stage audit](plus-source-stage-coverage-audit-2026-10-16.md).
