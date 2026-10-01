# N4 CLI synthetic fixture

`smooth_bias.nii.gz` is a 16×16×16 analytic scalar phantom. Its values are
`100 * exp(0.25*x + 0.16*y - 0.12*z)` on each axis sampled uniformly from
[-1, 1]. The affine has anisotropic 1.1×1.2×1.3 mm voxel sizes and a left-to-
right first voxel axis so that the CLI test also checks the required RAS
reorientation warning and output orientation.

The fixture exercises the command/API boundary only. N4 is supplied by
ANTsPy; the CLI comparison uses a short explicit three-iteration configuration
to keep the subprocess check quick and deterministic. The separate N4 API test
checks bias reduction using the selected wrapper defaults.
