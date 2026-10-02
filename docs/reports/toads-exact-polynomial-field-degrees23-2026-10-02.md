# TOADS exact polynomial field degrees 2 and 3 (P18cyau)

## Ported behavior

`fit_toads_exact_chebyshev_field` now supports the selected global exact-membership
Chebyshev correction fit at degrees 1, 2, and 3 (4, 10, and 20 terms). The
implementation preserves the selected `q=2` membership weights, the
constructor's IMAGE-to-OPTIMIZED_IMAGE switch for the optimized tissue atlas,
class-zero exclusion, per-category `PowerTable` weighting, every-third-voxel
sampling, source-order float32 terms, full-volume scaling, full-grid field
evaluation, and mask-mean shift. The prior degree-one public function remains
as a backward-compatible entry point.

## Synthetic source comparison

Committed test data in
[`tests/data/synthetic/toads_exact_polynomial_field_degrees23/`](../../tests/data/synthetic/toads_exact_polynomial_field_degrees23/)
is a 21³ RAS folded-cortex shell with a smooth positive image, eleven
normalized membership maps, and a unity starting field. Image payloads are
stored as `.nii.gz`. The Java formula harness uses the preserved Jama and
PowerTable reference sources and the selected method's operation order; this
is a method-formula comparison, not execution of the complete VandyCRUISE
wrapper.

The single-channel degree-2 and degree-3 normal matrices and right-hand sides
match the harness exactly. Degree-3 coefficients differ by at most `2.33e-11`
due to double-precision solver rounding. Complete single-channel fields
compare within `2e-6` absolute error, remain exactly one outside the correction
mask, and have within-mask float64 means within `2e-6` of one. A follow-on
two-channel check applies degrees 2 and 3 in source channel order and compares
each channel's normal system and output field against its own Java reference.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_exact_polynomial_field_degrees23.py tests/test_toads_multichannel_exact_field_refit.py tests/test_toads_post_thinning_exact_field.py tests/test_toads_post_thinning_exact_iteration.py tests/test_toads_p18cza_carried_degree_iteration.py
```

Result: **13 passed, 0 failed, 0 skipped in 77.60 seconds**. Implementation
and fixture commit: `4b1895e6300ecb9a597bf3ef2a164fabef0f0a5a`.
The Java harness SHA-256 is
`3dcca989aa148e126dd6a53101fd0f6d5f59a23fbc44686e97dd8f0f1283106f`.

The full repository suite has not yet been rerun on this commit, so P18cyau
remains Started and this is not full TOADS or MaCRUISE workflow evidence.
