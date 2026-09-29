# Full test suite: LesionTOADS smoothing and prior membership integration

- Date: 2026-09-29
- Tested implementation commit: `eb7e045a8927e3394a49e2f08c9352eae1f5cde0`
- Command: `.venv/bin/pytest -q`
- Result: **452 passed in 57.47 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes two selected-Java one-voxel compositions of smoothing,
shape/similarity prior updates, direct/prior membership inversion, lesion or
black-hole distance modulation, denominator normalization, and local energy.
The cases cover LESION/direct with ventricle and brainstem factors, and
BLACKHOLE/prior with the additional black-hole lesion-class prior and
GM-zero/brainstem modulation. All neighbor and shape/similarity arrays are
synthetic inputs; this is not a whole-volume segmentation validation.
