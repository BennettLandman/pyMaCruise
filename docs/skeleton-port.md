# ThinACE port progress

The selected source is Blake Lucas's VandyCRUISE `ThinACE.java`, called by the
MaCRUISE ACE adapter after ACE membership cleanup. Attribution and source
notices are retained in `licenses/upstream/` and `license_info.md`.

## Ported behavior

`src/pymacruise/skeleton.py` implements the source's thinning sequence:

1. Convert positive input values to object membership and clear the volume
   border. Build a signed binary level set and source-style distance priorities
   capped at four voxels.
2. Initialize the interior 6-neighbor surface list.
3. Classify points using a 3x3x3 local foreground 18-connectivity and
   background 6-connectivity test. Preserve source `surface` points and defer
   non-simple points to later waves.
4. Remove simple points in ascending distance order, retaining the source's
   endpoint rules based on 26-neighbor counts and eight octant checks. Add
   newly exposed 6-neighbors for subsequent waves.
5. Encode remaining skeleton membership as 254.

The public array API is `thin_ace`; `thin_ace_image` provides NIfTI input,
non-RAS warning/reorientation, and RAS output.

## Synthetic evidence

The source comparison harness compiles the original selected `ThinACE.java`
with small test-only JIST data-structure stubs; it runs the real algorithm
without a MIPAV/JIST application install. See
[`tests/java_reference/README.md`](../tests/java_reference/README.md) and
[`tests/data/synthetic/skeleton/README.md`](../tests/data/synthetic/skeleton/README.md).
Committed analytic `.nii.gz` inputs and Java outputs cover a line, cube,
nested GM/WM ellipsoid, and symmetric spheres of radii 7, 8, 9, 10, 11, 12,
14, and 16 voxels. Python output is voxel-exact against Java for every case
(Dice 1.0); symmetric mesh Hausdorff and mean nearest-surface distances are
both 0 mm for these equal masks. The sphere cases exercise equal-distance
priority ties at several sizes.

The Python port reproduces JIST heap insertion/removal tie behavior and the
source connectivity routines' neighborhood seed domains. Its distance field
uses the bundled MaCRUISE `DistanceFieldyh.java` reference, whose upstream
header identifies it as identical to the `DistanceField` class called by
ThinACE. NIfTI adapter tests confirm RAS output; a non-RAS input warns and is
reoriented. The project's acceptance limits are Dice >=0.97, symmetric
Hausdorff <=1 mm, and symmetric mean nearest-surface distance <=0.01 mm; the
committed source-reference fixtures exceed these limits exactly. Regenerate
goldens with `python scripts/generate_java_thinace_reference.py` (requires a
JDK and the project's Python dependencies).
