# Full test suite — AtlasEM initialization and output slices

- **Tested commit:** `1317cea7` (`Verify source-faithful AtlasEM zero prior behavior`)
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-16-atlasem-init-2 .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,606 passed, 4 skipped, 0 failed; 18 expected orientation warnings.
- **Elapsed:** 745.74 s (12 min 25.74 s).
- **Peak sampled RSS:** 14.372 GB; the 16 GB watchdog did not trigger.
- **Selected TOADS source-method estimate:** unchanged at 3,511–4,911 of 13,651 Java lines (26–36%). This checkpoint added AtlasEM initialization/membership/output slices and their source harnesses; it did not add TOADS methods. The estimate is source-method representation, not runtime coverage or end-to-end workflow completion.

The new AtlasEM methods use deterministic synthetic arrays and dependency-free
Java excerpts pinned to the supplied source hashes. This run does not validate
the missing `ibsr8obj` AtlasEM resource bundle, atlas registration, complete
AtlasEM evolution, pre-LesionTOADS Dura evolution, or whole PLUS equivalence.
See [P44b source-stage audit](plus-source-stage-coverage-audit-2026-10-16.md).
