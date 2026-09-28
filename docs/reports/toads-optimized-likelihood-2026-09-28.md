# Optimized TOADS tissue likelihood comparison

**Date:** 2026-09-28  
**Roadmap:** P18ae complete  
**Validated implementation/fixture commit:** `92550a4eeb11b2b049823a60b022581b0db0c84b`

## Source slice

`src/pymacruise/toads_math.py::toads_optimized_tissue_data_term` implements
the selected `MultipleToadSegmentation.computeExactMemberships` voxelwise
intensity term for the optimized CSF, GM, and WM class models. It applies the
optional multiplicative intensity field, PSF bias, per-modality residual and
squared distance, source-compatible `PowerTable` lookup, modality weight,
class variance weight, and image variance weight. Arithmetic stages that are
float in the source are explicitly retained as float32; the lookup helper
reuses NumPy power after reproducing the Java table index.

The comparison harness uses the committed actual `PowerTable.java` class and a
source-formula excerpt from `MultipleToadSegmentation`. SHA-256 values:

- `MultipleToadSegmentation.java`: `739a770b866dcc76b8198d795a06f2ae90eecaa437bed50c2fe151bf72bf1fb1`
- `PowerTable.java`: `2c75d78ff33efd11531a3227399acc6fb2a421e2b4205ad6c77e6533eadd50ae`

## Synthetic fixture and result

The committed `tests/data/synthetic/toads_optimized_likelihood/` fixture has
three 7×6×5×3 float32 RAS `.nii.gz` arrays: normalized input channels, a
multiplicative field, and PSF bias. The affine is 1 mm isotropic RAS. The
generator recipe and expected-value checks are deterministic and contain no
human imaging data.

The Java harness emits 630 values (210 voxels × three optimized tissue
classes) for each of two modes, with the field disabled and enabled. Both
include PSF bias. The focused test checks every voxel/class output in both
modes with `rtol=0`, `atol=1e-7`; all passed. The field-on and field-off results
also differ, confirming that the field is active in the calculation.

## Scope limits

This is the per-voxel tissue intensity data term only. It does not claim a
complete `computeExactMemberships` port: smoothing, atlas priors, outlier or
lesion mixtures, membership normalization, update ordering, topology gates,
and deformable TOADS evolution remain open. Reference results are formula
harness outputs, not full Java TOADS segmentation equivalence.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **271 passed in 59.58 seconds**
- Pytest: 9.1.1; Python: 3.13.2
- Tested implementation commit: `92550a4eeb11b2b049823a60b022581b0db0c84b`
- Run date: 2026-09-28
