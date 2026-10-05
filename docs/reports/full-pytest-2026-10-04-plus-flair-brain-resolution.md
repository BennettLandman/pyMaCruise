# Full repository pytest — brain-resolution plus FLAIR checkpoint

- Tested tree: `826e6ae65c3515c782c8b4788ad3a9f726a0c543` (P41d implementation/test commit `9e0654f817948b95cd5966bb436816139d0c1d6b`).
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,168 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 877.65s (14:37)**.
- Peak sampled RSS: **13.971 GB**; the 24 GB watchdog did not terminate pytest.

This run includes the committed brain-resolution synthetic FLAIR-to-T1 checkpoint test. It validates a finite RAS warped image on the fixed T1 grid, persisted transforms and manifest, and improvement in synthetic whole-grid correlation. That metric is a smoke/geometry check; it is not physical transform-accuracy acceptance or NiftyReg/Aladin numerical equivalence. No human MRI was used.
