# Full test suite: LesionTOADS black-hole distance modulation

- Date: 2026-09-29
- Tested implementation commit: `c5493e84aa72d03bd951194aad0a91eec0ee3e39`
- Command: `.venv/bin/pytest -q`
- Result: **442 passed in 81.37 seconds**
- Warnings: two expected warnings for synthetic Dura inputs that are reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes fourteen selected-Java formula cases for the
`BLACKHOLE` distance-modulation branch. Lesion membership,
`blackHoleMem`, and `lesionRelation` float32 output bits match. The cases
exercise the source's ventricle branch, which leaves `blackHoleMem`
unchanged; the GM branch, which scales both auxiliary maps; independent
brainstem modulation; zero-distance `spread`; `-1` sentinels; zero maxima;
and GM enable/disable. The cases use already-computed scalar distances; they
do not validate distance-map construction or normalized voxel memberships.
