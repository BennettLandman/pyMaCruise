# Full repository pytest — plus LesionTOADS image preparation

- Tested implementation/fix commit: `49a896b4fd6fe5622a2358c795da3d0bdc7ce4f8`.
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,171 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 870.25s (14:30)**.
- Peak sampled RSS: **14.424 GB**; the 24 GB watchdog did not terminate pytest.

The run includes P41e: source-order rigid FLAIR-to-T1 registration, N4 correction
of T1 and registered FLAIR, application of the required supplied brain mask,
and restartable RAS modality outputs. A first run at implementation commit
`5f67554` exposed a local `N4Config` import-shadowing regression in unrelated
CLI paths; a focused regression group passed after the fix before this full
rerun. The selected TOADS source-line estimate did not change in this plus
workflow cycle. No human MRI was used.
