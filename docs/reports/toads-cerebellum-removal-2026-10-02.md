# Selected second-pass cerebellum removal (P18cybr)

## Ported behavior

The MATLAB MaCRUISE wrapper applies cerebellum removal independently to the
second-pass TOADS GM, WM-fill, and WM-mask maps. It transfers MaSeg labels to
each map with `ToadsCruiseSpaceToOriginalSpace` and zeros labels 35, 38–41,
and 71–73. The Python API uses the existing nearest-neighbor physical-grid
transfer helper, preserves each tissue map's dtype, and returns a changed-voxel
count. The CLI writes all three cleaned maps as RAS `.nii.gz` and can record
input, parameters, and outputs in a provenance JSON file.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_cerebellum.py
```

Result: **2 passed, 0 failed, 0 skipped in 3.73 seconds**.

The committed 9×8×7 RAS phantom contains each removed label once, a retained
label, and background. API checks confirm exactly eight selected voxels are
zeroed while all other values, image geometry, and map dtype are preserved.
The CLI test writes three native-RAS float32 products and provenance.

The source JIST transfer uses orientation-based `ReorientVolume.solve` with
unchanged resolution. The Python transfer uses physical affines; this test
establishes exact behavior for identical RAS grids, not arbitrary differing
origins or spacing. Full pipeline execution is not claimed. P18cybr remains
Started.
