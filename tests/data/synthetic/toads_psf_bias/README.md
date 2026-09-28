# TOADS PSF ridge-bias phantom

This deterministic synthetic fixture contains 9×8×7×2 input and field arrays
and a 9×8×7 integer segmentation, all stored as RAS `.nii.gz` with a 1 mm
isotropic affine. `scripts/generate_toads_psf_fixture.py` defines the voxel
values and labels. It includes a planar tissue boundary and a small island so
the 3×3×3 boundary condition and center-channel calculation are exercised.

`java_reference.json` contains the source-formula PSF bias on the interior
grid, with and without multiplication by the intensity field. The excerpt
matches `MultipleToadSegmentation.isBoundary` and the PSF precomputation in
`computeExactMemberships`; the source file SHA-256 is recorded in the JSON.
Java float32 operation order is retained, including small roundoff residuals
from summing the center and all 26 surrounding stencil contributions.

Regenerate NIfTI inputs with `scripts/generate_toads_psf_fixture.py` and the
Java reference with `scripts/generate_java_toads_psf_reference.sh`.
