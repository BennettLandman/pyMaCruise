# Selected second-pass GM-derived TOADS crop (P18cybs)

## Ported behavior

The selected base pipeline derives one crop box from positive second-pass GM
support, expands the inclusive bounds by the configured edge width (5 voxels
by default), then reuses those exact crop parameters for WM-fill and WM-mask.
The Python API and `crop-tissue-maps` CLI use the existing `CubicVolumeCropper`
adapter, retain scalar dtypes, pad outside-grid voxels with zero, and preserve
the crop affine in RAS.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_tissue_crop.py tests/test_cropping.py tests/test_plus_lesion_transfer_crop.py
```

Result: **7 passed, 0 failed, 0 skipped in 33.37 seconds**.

The committed asymmetric 13×12×11 RAS fixture has a folded/holed GM support
and coordinate-coded target maps. The Python stage checks expected bounds,
dimensions, target voxel order, and RAS crop transforms. The related
plus-lesion integration test compiles and runs the selected identity-RAS Java
reference harness over eight companion arrays, including the same shared
cropper operation.

This verifies the selected crop on identical RAS grids. JIST transfer to an
arbitrary SPECTRE/reference grid remains a separate limitation. Topology
correction, ACE, GVF, and TGDM after this crop are not included in this stage.
P18cybs remains Started.
