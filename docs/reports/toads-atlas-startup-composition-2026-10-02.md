# P18cyad: TOADS atlas assets through initial registration

## Scope

Implementation commit: `e80ef81a40195b195eae32c31a4033b8e95262c3`
(`Compose TOADS atlas assets with initial registration`).

`register_toads_from_atlas_assets` wires an atlas descriptor or parsed atlas,
RAS topology and shape maps, one modality's explicit descriptor intensity
priors, the one-channel T1/user-mask startup, descriptor-order shape stacking,
and the existing exact-membership rigid pyramid. Pyramid depth, iteration cap,
and modality remain explicit parameters.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_atlas_startup.py tests/test_toads_alignment_startup.py tests/test_atlas_assets.py
```

Result on the exact implementation commit: **13 passed, 0 failed, 0 skipped in
1.24 seconds**. The new compact synthetic two-class descriptor and compressed
NIfTI atlas assets verify prior parsing, class ordering, RAS source-orientation
metadata, initial membership class count, one scale-2 registration level, and a
finite resulting transform. No real imaging data was read.

## Limits

This is adapter composition only. It does not independently compare the full
startup/registration chain against Java or reproduce every state transition in
`initShapeRegistration` and `updateShapeRegistration`; it does not perform
TOADS iteration or claim workflow equivalence. The full test suite was not run
for this localized change.
