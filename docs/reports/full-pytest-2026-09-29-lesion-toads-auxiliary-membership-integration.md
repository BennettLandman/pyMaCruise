# Full test suite: LesionTOADS auxiliary-membership integration

- Date: 2026-09-29
- Tested implementation commit: `edad578985013fe3a747a4fad42df11c22a73013`
- Command: `.venv/bin/pytest -q`
- Result: **455 passed in 56.82 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes two selected-Java one-voxel compositions that place a
LESION or BLACKHOLE auxiliary map and an OUTMASK/outlier map in the same
normalization denominator. Class/auxiliary memberships, relation,
denominator, and local energy match. Tests include GM-zero/brainstem,
ventricle precedence with the source's unchanged black-hole map, and
direct/prior OUTMASK mixing. Input memberships and energies are preconverted.
