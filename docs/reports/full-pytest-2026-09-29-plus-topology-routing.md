# Full test suite: MaCRUISEplus MaWM topology routing

- Date: 2026-09-29
- Tested implementation commit: `68f013ffb24d8f8e7f65c9124c9627cb76b4968a`
- Command: `.venv/bin/pytest -q`
- Result: **473 passed in 56.75 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The plus topology adapter receives the committed 9³ RAS topology scalar field,
paint mask, and MaSeg-derived WM image. With the selected refinement flag, it
routes MaWM into `topology_correct(..., multi_atlas_wm=...)` and matches the
direct result, including the refined RAS output. It fails clearly if MaWM is
missing while refinement is enabled. With refinement disabled, it omits MaWM,
returns no atlas-refined product, and matches the base topology output.

Focused plus/topology tests: 16 passed. The three `.nii.gz` inputs, generator,
and hashes are committed in `tests/data/synthetic/plus_topology/`. This
validates option and input routing to the source-checked refinement; it does
not complete plus workflow orchestration or full cortical topology outputs.
