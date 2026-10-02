# TOADS initial image-variance composition (P18cyan)

P18cyan adds the selected `computeImageVariance()` startup operation to the
atlas-backed registration-to-segmentation initialization chain. After exact
memberships and assigned-class deviations are available, the adapter scans the
full normalized image grid for optimized CSF/GM/WM classes, excludes MASK, and
computes normalized inverse-variance channel weights and image-level standard
deviations. Startup correction fields are unity, so the source-equivalent
calculation receives the normalized images directly.

The component kernel retains source float32 accumulation, x/y/z traversal,
membership-squared weighting with `1e-30f` stabilization, fourth-power
residuals, and the source zero-residual INF/NaN behavior. This composition does
not yet calculate centroid similarities or run the exact membership/energy
sweep.

Focused validation on synthetic data:

- Command: `./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_initial_relation_state.py tests/test_toads_relation_list.py tests/test_toads_boundary_distances.py tests/test_toads_closest_labels.py tests/test_toads_relation_state_refresh.py tests/test_toads_optimized_class_variance.py tests/test_toads_optimized_variance.py tests/test_toads_object_image_variance.py`
- Result: **45 passed, 0 failed, 0 skipped in 8.01 seconds**.
- Implementation commit: `5a4aca2567bcfa138020ab3d410337e1b55cee20`.

The latest full-suite result remains the separate checkpoint on
`ee350877e61f97adb050d185bfb34f7bfd067ade`; this focused validation does not
replace it or claim full TOADS workflow equivalence.
