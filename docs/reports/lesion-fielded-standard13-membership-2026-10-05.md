# Fielded standard-13 LesionTOADS membership update — 2026-10-05

## Implementation

The selected Java loop can attach one multiplicative correction field to each
image channel. Its membership calculation uses `field[channel] * image[channel]`
before active OBJECT and LESION residual calculations; MASK mismatch tests the
corrected channel value against zero. The Python standard-13 voxel and volume
passes now accept and validate those fields, and the next GM-enabled membership
adapter and reusable descriptor loop forward them. Existing callers without
fields retain their prior behavior. The optimized-CSF branch continues to use
the selected source's active scalar outlier fallback.

This implements application of already-attached fields during a membership
pass. It does not refit fields at each loop iteration. The adaptive degree gate
and selected global Chebyshev refit exist as separate tested components, but
their integration into repeated LesionTOADS execution is still open.

Implementation commit: `65ed0fe5641abf7a3b2ba2cf8344e2ccc83e334f`.

## Evidence

Focused component/formula command:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q --tb=short \
  tests/test_lesion_standard13_voxel.py \
  tests/test_lesion_standard13_volume.py \
  tests/test_lesion_post_thinning_gm_membership.py \
  tests/test_lesion_descriptor_iteration_slice.py \
  tests/test_lesion_toads_data_terms.py \
  tests/test_lesion_active_class_energy.py \
  tests/test_lesion_toads_field_data_terms.py \
  tests/test_lesion_toads_fielded_data_energy_membership.py
```

Result: **21 passed, 0 failed, 0 skipped in 11.97s**. This includes the
selected Java-backed fielded data-term references and new identity-field,
non-unit-field, and validation checks at voxel and volume scope.

Registered-loop integration command:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q --tb=short \
  tests/test_lesion_registered_repeated_iteration.py::test_registered_loop_driver_matches_two_manual_updates
```

Result: **1 passed, 0 failed, 0 skipped in 16.08s**. A compact RAS synthetic
case carries non-unit T1/FLAIR fields through the next statistics and
membership transition; the resulting objective differs from the uncorrected
control. Full suite on implementation commit `2b90545a52ace6a3d518925a380c9c34ccac88a0`: **1,200 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,239.07s**, peak sampled RSS **12.993 GB** under the 24 GB watchdog; see [full-suite report](full-pytest-2026-10-05-lesion-adaptive-fields.md).

## Limits

The fixture is synthetic. It validates Python composition around separately
Java-checked data-term behavior; it is not a full Java `computeMemberships`
comparison or real-image result. Adaptive field fitting, persistent field
state at restart boundaries, source-composite equivalence, and human MRI
validation remain open.
