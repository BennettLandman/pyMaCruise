# Full pytest run after JIST-compatible Dura transfer mode — 2026-10-11

- Tested implementation commit: `eb10ccfe` (`Add opt-in JIST-compatible Dura transfer mode`)
- Command: `PYTHONPATH=src PYTHONNOUSERSITE=1 MPLCONFIGDIR=build/mpl-cache .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- .venv/bin/pytest -q`
- Result: **1,559 passed, 4 skipped, 0 failed**, 11 expected RAS reorientation warnings
- Runtime: **812.06 s** (13 min 32 s)
- Peak sampled RSS: **11.591 GB**; 16 GB watchdog did not terminate the run

The suite includes default-on optional RAS regridding, physical-affine Dura
transfer, the opt-in `jist-unchanged` transfer/crop path, and its synthetic
full raw-input `plus-full-run` integration case. The JIST compatibility route
is reproducible and restartable; the P41z downstream anatomical quality audit
still fails the agreed 97% support gate and is reported separately in the
[transfer downstream report](p41z-jist-transfer-downstream-2026-10-11.md).

This is synthetic and unit/integration validation. It is not full Java
MaCRUISEplus equivalence or real-MRI acceptance.
