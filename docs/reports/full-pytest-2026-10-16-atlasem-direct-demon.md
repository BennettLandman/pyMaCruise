# Full test suite — AtlasEM direct Demon optimizer

- **Tested code commit:** `e11753d350adc9db80a2e8619b9280f7f83bb4fc`
- **Command:**
  `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLCONFIGDIR=build/full-suite-mpl-2026-10-16-atlasem-demon-java .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 2 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- **Result:** 1,648 passed, 4 skipped, 0 failed; 21 expected orientation warnings.
- **Elapsed:** 736.03 s (12 min 16.03 s).
- **Peak sampled RSS:** 14.206 GB; the 16 GB watchdog did not trigger.

The suite includes direct selected-Java `DemonToadDeformableAtlas` rigid
registration comparisons for a synthetic one-step case and a three-step cap.
The six final transform parameters match bit for bit. This whole-suite result
does not supply the missing selected `ibsr8obj` descriptor/payload or validate
complete AtlasEM/PLUS behavior.
