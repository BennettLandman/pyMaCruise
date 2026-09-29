# Full test suite: LesionTOADS six-face smoothing energy

- Date: 2026-09-29
- Tested implementation commit: `bba409dc9c892262067f8ceaeca06c5582d58e7b`
- Command: `.venv/bin/pytest -q`
- Result: **449 passed in 52.76 seconds**
- Warnings: two expected warnings for synthetic Dura inputs that are reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes four selected-Java one-voxel smoothing cases for
`OBJECT`, `LESION`, `BLACKHOLE`, and `OUTMASK`. Full membership data-term
increments, half-weighted local energy increments, and intermediate neighbor
square sums match. The cases verify face/class traversal order,
current-class exclusion, lesion/black-hole cross-coupling, and the selected
source's self-neighbor ordering. This does not validate volume-wide smoothing
or segmentation evolution.
