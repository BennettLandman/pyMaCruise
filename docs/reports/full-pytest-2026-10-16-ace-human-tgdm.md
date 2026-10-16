# Full test suite — ACE wrapper and human BASE TGDM evidence checkpoint

- **Tested tree commit:** `857370c002164a5775633117059fe35ecdd7aaa8`
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-16 .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,582 passed, 4 skipped, 0 failed; 18 expected orientation warnings.
- **Elapsed:** 748.41 s (12 min 28 s).
- **Peak sampled RSS:** 7.205 GB; 16 GB watchdog did not trigger.
- **Source-method estimate:** unchanged at 3,511–4,911 of 13,651 selected TOADS Java lines (26–36%). This checkpoint added tests/evidence only; it did not change source-method mappings, so test count is not used as a port-coverage proxy.

The suite ran after the selected human BASE attempt-02 TGDM Java comparison and
documentation reconciliation. It includes the source-faithful ACE float-call
argument correction and the current human BASE workflows. It does not replace
the selected-stage Java comparisons or establish whole-workflow equivalence.
