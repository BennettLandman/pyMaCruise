# Full test suite: LesionTOADS data-energy membership integration

- Date: 2026-09-29
- Tested implementation commit: `990407f49bf7ca271567f032a0d7d62a0bddfb06`
- Command: `.venv/bin/pytest -q`
- Result: **456 passed in 64.21 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes two selected-Java one-voxel compositions beginning with
two-channel unfielded image residuals for LESION and BLACKHOLE classes. Data
energies match before smoothing, shape/similarity priors, inverse conversion,
distance modulation, normalization, and local energy. Spatial stencils and
shape/similarity vectors are fixed synthetic inputs; there is no volume run.
