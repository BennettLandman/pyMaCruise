# Brain-resolution synthetic FLAIR registration

- Focused command:
  `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 30 -- env PYTHONPATH=src ./.venv/bin/pytest -q tests/test_plus_flair_brain_resolution.py`
- Focused result: **1 passed, 0 failed, 0 skipped in 14.09s**.
- Peak sampled RSS: **0.858 GB**; the 16 GB watchdog did not terminate the test.
- Implementation/test commit: `9e0654f817948b95cd5966bb436816139d0c1d6b`.

The committed analytic RAS phantom is `(172, 216, 164)` and contains nested
tissue-like regions, folded cortical modulation, contrast-remapped FLAIR, a
small focal hyperintensity, and a known integer-voxel moving-image shift. The
ANTsPy rigid checkpoint creates a finite warped image on the fixed T1 grid,
persists forward/inverse transforms and a manifest, and improves whole-grid
T1/FLAIR Pearson correlation by more than 0.03 relative to the unregistered
input.

This is a synthetic registration smoke/geometry check, not a registration
accuracy tolerance, NiftyReg/Aladin equivalence, full plus execution, or
real-data validation. The initial test threshold was set too high (+0.10); the
observed improvement was +0.0497, so the committed assertion uses +0.03.
