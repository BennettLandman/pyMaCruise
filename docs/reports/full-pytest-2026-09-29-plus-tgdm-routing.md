# Full test suite: MaCRUISEplus cropped TGDM routing

- Date: 2026-09-29
- Tested implementation commit: `2e9503f70cb8fe9f2169d6cd065eaac5f67c932e`
- Command: `.venv/bin/pytest -q`
- Result: **480 passed in 79.51 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

`run_plus_tgdm_stage` routes the cropped corrected level set, updated WM and
GM, MaACE-enhanced GM, thinned MaACE skeleton, and plus GVF volume to the
existing nested TGDM image adapter. It obtains the selected WMIsoLevel through
`plus_tgdm_parameters` and preserves the cropped RAS grid. On the committed
25³ synthetic TGDM phantom with WMIsoLevel set to 0.43, all six image products
match a direct `tgdm_nested_images` call with the same mapped parameters;
misaligned GVF input is rejected.

Focused TGDM tests: 22 passed in 37.52 seconds (plus adapter and TGDM image
tests). This validates plus input and parameter routing into the existing
surface engine. Inverse cropping, output persistence, and full plus workflow
orchestration remain separate work.
