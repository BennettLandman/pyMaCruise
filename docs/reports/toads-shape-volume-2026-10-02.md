# TOADS transformed atlas shape volume (P18cyas)

P18cyas adds the selected atlas-shape sampling stage over the current
post-registration RAS crop. It uses descriptor class order, the final rigid
transform and class bounds from the registration adapter, signed subject
spacing, descriptor template spacing, source linear-closest interpolation,
shape-slope recovery (`2.0`), and the class-zero outside-bound default. The
vectorized output is `(x,y,z,class)` and is exposed on
`ToadsAtlasAlignmentRun`.

The anisotropic synthetic test compares every vectorized voxel to the existing
selected pointwise `sample_toad_shape_prior` adapter, with absolute tolerance
`3e-7`. This covers the small float difference between regular point mapping
and the source's precomputed matrix route. It also checks class-specific
exclusive upper bounds, empty-support defaults, and integration with the
synthetic atlas startup. The transformed volume is now available from startup;
P18cyaq/P18cyar still accept transformed shapes as an explicit input and do
not yet consume this returned volume.

- Command: `./.venv/bin/pytest -q tests/test_toads_shape_volume.py tests/test_toads_shape_sampling.py tests/test_toads_atlas_startup.py tests/test_toads_post_registration_template_refresh.py tests/test_toads_registration.py`
- Result: **35 passed, 0 failed, 0 skipped in 7.77 seconds**.
- Implementation commit: `8088396d727e38187bc8bddfdd6db08d01198dda`.
