# Full test suite — AtlasEM to pre-LesionTOADS Dura handoff

- **Tested commit:** `1ef97abc` (`Connect AtlasEM outputs to pre-LesionTOADS Dura`)
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-16-atlasem-dura-handoff .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,619 passed, 4 skipped, 0 failed; 18 expected orientation warnings.
- **Elapsed:** 779.66 s (12 min 59.66 s).
- **Peak sampled RSS:** 5.737 GB; the 16 GB watchdog did not trigger.
- **Selected TOADS source-method estimate:** unchanged at 3,511–4,911 of 13,651 Java lines (26–36%). The current work adds AtlasEM/Dura slices but no methods to the selected TOADS estimate.

The added integration case runs source-matched AtlasEM output formation into a
synthetic 17³ RAS pre-LesionTOADS Dura pass. This is synthetic component
composition, not selected-atlas, original-workflow, or real PLUS validation.
The `ibsr8obj` source resource bundle and atlas-registration stage remain
unavailable. See [P44b source-stage audit](plus-source-stage-coverage-audit-2026-10-16.md).
