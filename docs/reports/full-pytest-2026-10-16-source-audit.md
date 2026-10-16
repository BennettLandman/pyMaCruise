# Full test suite — source-reference packet and PLUS source audit

- **Tested tree commit:** `ca25a55661e863959e5b1ce2220815e7abe4f5b0`
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-16-source-audit .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,591 passed, 4 skipped, 0 failed; 18 expected orientation warnings.
- **Elapsed:** 782.18 s (13 min 02 s).
- **Peak sampled RSS:** 11.164 GB; the 16 GB watchdog did not trigger.
- **Source-method estimate:** unchanged at 3,511–4,911 of 13,651 selected TOADS Java lines (26–36%). No TOADS methods were added or remapped in this checkpoint.

This full suite covers the committed source-output auditor, BASE/PLUS
selected-mask driver wiring guards, atlas-resource preflights, and source-stage
coverage audit. It does not parse or run MATLAB, Java plugins, or Aladin and
does not validate AtlasEM/Dura output parity. See the
[PLUS source-stage audit](plus-source-stage-coverage-audit-2026-10-16.md) and
[source-reference packet](../../source-reference-run/README.md) for those
remaining boundaries.
