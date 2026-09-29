# Full test suite: MaCRUISEplus MaACE stage wiring

- Date: 2026-09-29
- Tested implementation commit: `89514b8d4886d668733af2e9ec204a6b641a99c9`
- Command: `.venv/bin/pytest -q`
- Result: **470 passed in 57.48 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The plus-specific ACE adapter accepts typed plus settings and passes the
selected `use_updated_csf_in_ace=3` placement and `Merge` behavior to the
existing MaACE implementation. With the committed 24³ RAS synthetic phantom,
enhanced GM agrees with the Java after-Gaussian reference within MAE
`1e-6` / maximum absolute error `5e-6`; cleaned and thinned skeleton masks
match exactly. The `update_Enhance=0` selection is observed. An enabled update
value rejects explicitly because that later branch is outside this adapter.

Focused plus ACE/reference checks: 5 passed. This connects the typed profile
to the validated ACE stage; the full MaCRUISEplus workflow and separate
MaACE paper boundary prior remain open.
