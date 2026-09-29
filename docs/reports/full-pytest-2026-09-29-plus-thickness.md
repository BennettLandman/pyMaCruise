# Full test suite: MaCRUISEplus thickness routing

- Date: 2026-09-29
- Tested implementation commit: `2113eb058ab4e5fd124b8994f95efc7693166acc`
- Command: `.venv/bin/pytest -q`
- Result: **484 passed in 90.08 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

`run_plus_thickness_stage` routes the uncropped inner/outer TGDM volumes to the
existing level-set subtraction function and the corresponding meshes to the
existing closest-point thickness function. On the committed nested-sphere
fixture, the volume contains only 2 mm differences, with mean 2 mm and sample
standard deviation 0 across all crossing voxels. Surface statistics and
per-vertex thickness match direct calls to `closest_point_thickness`; outputs
are RAS+ millimetres.

Focused tests: 9 passed in 4.59 seconds, including the existing thickness and
surface suites. This validates routing into the existing source-mapped
thickness algorithms. Surface thickness remains subject to the documented
JIST-versus-Lewiner mesh substitution risk.
