# Full test suite: MaCRUISEplus WMIsoLevel routing

- Date: 2026-09-29
- Tested implementation commit: `cef0c45023013c3fd5ba0d51a9ceec2ecc0a8510`
- Command: `.venv/bin/pytest -q`
- Result: **468 passed in 57.98 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

`plus_tgdm_parameters` maps the selected MaCRUISEplus `WMIsoLevel=0.5` to
`TGDMBaselineParameters.white_matter_iso`. A second check sets a synthetic
nondefault isolevel and verifies that only that TGDM field changes; all other
selected defaults remain equal to the baseline object. Focused plus/TGDM
tests: 29 passed. This validates parameter routing, not a MaCRUISEplus volume
run or nested-surface output.
