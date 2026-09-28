# ThinACE synthetic reference fixtures

All images are synthetic, 1 mm isotropic RAS `.nii.gz` files. Inputs are
analytic binary objects except for `nested_gm_wm_input.nii.gz`, which encodes
an analytic nested gray/white matter ellipsoid using values 128 and 255.
Corresponding `*_thinace_java.nii.gz` files are outputs from the selected
VandyCRUISE `ThinACE.java` implementation in
`tests/java_reference/upstream/edu/jhu/ece/iacl/algorithms/ace/`.

Cases include a one-voxel line, solid cube, nested GM/WM ellipsoid, and
symmetric spheres with radii 7, 8, 9, 10, 11, 12, 14, and 16 voxels. The
sphere sweep checks source-priority ties as the object size grows. Each output
is compared voxel-for-voxel; `tests/test_skeleton_reference.py` also checks
Dice and symmetric mesh surface distances against project tolerances.

Regenerate fixtures and `SHA256SUMS` from the repository root with:

```sh
python scripts/generate_java_thinace_reference.py
```

This requires a JDK and the project Python dependencies. Java class files and
binary interchange files are transient under ignored `build/`; only the
compressed NIfTI images, metadata, and checksums are committed. Source and
license provenance is described in `tests/java_reference/README.md` and
`license_info.md`.
