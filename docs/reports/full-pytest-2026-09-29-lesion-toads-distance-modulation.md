# Full test suite: LesionTOADS distance modulation

- Date: 2026-09-29
- Tested implementation commit: `c93ec3978753acf857a656c0d4792228d176c90d`
- Command: `.venv/bin/pytest -q`
- Result: **440 passed in 69.12 seconds**
- Warnings: two expected warnings for synthetic Dura inputs that are reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes the P39z selected-Java formula comparison. Twelve scalar
cases match lesion membership and `lesionRelation` float32 output bits across
ventricle/GM precedence, zero-distance `spread`, absent-distance sentinels,
disabled maxima, GM enable/disable, and independent brainstem modulation.
The selected cases exercise already-computed scalar distances; they do not
validate construction of the distance maps, black-hole distance modulation,
or normalized voxelwise memberships.
