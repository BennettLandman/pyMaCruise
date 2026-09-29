# Full test suite: LesionTOADS smoothing-to-membership integration

- Date: 2026-09-29
- Tested implementation commit: `eb0b7cc5d13a959608ad7acabc63675a494874cd`
- Command: `.venv/bin/pytest -q`
- Result: **450 passed in 57.57 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes two selected-Java one-voxel compositions from the six-face
smoothing update through class/auxiliary data-term inversion, relation
scaling, lesion or black-hole distance modulation, denominator
normalization, and local energy. The cases cover LESION/direct and
BLACKHOLE/prior branches, ventricle/GM selection, zero GM distance, and
brainstem distance. Atlas-prior energy construction and volume-wide updates
are not included.
