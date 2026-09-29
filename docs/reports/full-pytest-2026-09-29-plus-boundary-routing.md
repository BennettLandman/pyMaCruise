# Full test suite: MaCRUISEplus boundary routing

- Date: 2026-09-29
- Tested implementation commit: `7741bb27db7be5a3d51c9adda9c510fb0c1d452c`
- Command: `.venv/bin/pytest -q`
- Result: **475 passed in 57.23 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

The plus boundary adapter applies the configured overlooked labels during the
selected 3D multi-atlas boundary construction, then clears boundary voxels
where the aligned TOADS WM mask is strictly greater than 0.01. The committed
13³ anisotropic RAS fixture includes label 120 and partial WM coverage. Tests
compare both adapter products to direct calls of the existing helpers, verify
that overlooking label 120 changes the raw boundary, verify WM exclusion, and
reject an unsupported erosion method.

Focused boundary tests: 13 passed. The segmentation and TOADS WM NIfTI inputs,
generator, provenance note, and SHA-256 manifest are committed in
`tests/data/synthetic/plus_boundary/`. This validates plus-stage routing and
composition of existing source-checked operations; it does not complete the
plus workflow or validate full cortical surfaces.
