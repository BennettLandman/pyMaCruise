# Full test suite: LesionTOADS membership integration

- Date: 2026-09-29
- Tested implementation commit: `a55cc41d6c451f8d3698d1e80ee5128c9a909ee1`
- Command: `.venv/bin/pytest -q`
- Result: **448 passed in 61.93 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes five selected-Java one-voxel formula compositions. Class
memberships, lesion/black-hole auxiliary memberships, `lesionRelation`,
denominator, and local energy match through inverse conversion, direct/prior
scaling, distance modulation, and final normalization. The cases cover mixed
class types, ventricle/GM/brainstem distance paths, and a zero denominator.
This does not validate smoothing, prior-energy construction, or a full
multi-voxel membership update.
