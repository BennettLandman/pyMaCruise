# Full test suite: LesionTOADS atlas-prior energy

- Date: 2026-09-29
- Tested implementation commit: `0a7fc6ca18fcc433568711839515f3e506201c10`
- Command: `.venv/bin/pytest -q`
- Result: **451 passed in 57.63 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes six selected-Java one-voxel cases for the
`priorCoefficient` source excerpt. Class, lesion, and black-hole prior sums,
full membership data-term additions, and half-weighted energy additions match
at the float32 bit level. The cases cover ordinary, lesion, black-hole,
optimized-lesion, and outmask types; optional lesion-class terms; a
black-hole lesion class; current-class exclusion; zero priors/shapes; and
float32-sensitive values. Shape and similarity arrays are supplied inputs.
