# Full test suite: LesionTOADS voxel membership normalization

- Date: 2026-09-29
- Tested implementation commit: `4d0a25af19d146823402bea1febd5dfc60868695`
- Command: `.venv/bin/pytest -q`
- Result: **447 passed in 67.49 seconds**
- Warnings: two expected warnings for synthetic Dura inputs that are reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes eight selected-Java one-voxel cases for the final
`computeMemberships` normalization excerpt. Denominator, class memberships,
lesion/black-hole/outlier auxiliary memberships, and local energy match at
the float32 bit level. Cases cover ordinary, lesion, black-hole,
outlier/outmask, nonpositive denominator, repeated auxiliary division,
optimized lesion, and float32-sensitive values. Inputs are pre-assigned
memberships and energies; the suite does not claim full voxelwise membership
estimation.
