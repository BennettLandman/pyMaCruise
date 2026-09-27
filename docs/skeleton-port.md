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

## Synthetic evidence and open checks

Tests check exact retention of a one-voxel line, volume reduction of a solid
cube, 18-connected output, border clearing, and RAS warning/output behavior.
The committed ACE concentric phantom inputs are under
`tests/data/synthetic/ace/`; `make_phantom.py` regenerates the `.nii.gz` images
and their checksums. The full offline test suite passes 85 tests.

The implementation is not yet marked validated: Java heap tie ordering has not
been reproduced byte-for-byte, the historical `DistanceField` source used by
ThinACE is not present as a Java file in the supplied checkout, and no trusted
Java result mask has yet been generated. Compare line, cube, sphere, and nested
GM/WM cases to the source runtime before closing P27; retain the project's Dice
and surface-distance thresholds.
