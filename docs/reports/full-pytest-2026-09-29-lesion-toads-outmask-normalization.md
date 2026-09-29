# Full test suite: LesionTOADS OUTMASK normalization integration

- Date: 2026-09-29
- Tested implementation commit: `899493c1c3999871b5c57f88a5df9a033f8701d9`
- Command: `.venv/bin/pytest -q`
- Result: **454 passed in 56.88 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The run includes three selected-Java one-voxel compositions of OUTMASK
mixture, outlier-map update, denominator accumulation, membership
normalization, and local energy. Direct/prior relations, zero/unit outlier
weights, and repeated OUTMASK class updates are covered. The fixture begins
with preconverted class/outmask values and does not validate their data terms.
