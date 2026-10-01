# TOADS inhomogeneity polynomial-basis fixture

`grid.nii.gz` is a synthetic 5 x 4 x 3 RAS grid with anisotropic voxel sizes
0.8 x 1.1 x 1.7 mm. `generate_fixture.py` recreates it. The Java harness reads
the NIfTI dimensions and evaluates the selected formula excerpt from
`InhomogeneityCorrection2.compute3DChebyshevCoefficients` at every voxel for
degrees 1, 2, and 3. Its source hash is pinned in `java_reference.json`.

All coefficients match the Java float32 results bit-for-bit. The tested basis
only is a reusable sub-operation; it does not fit polynomial coefficients,
solve the weighted least-squares system, produce/apply the multiplicative
field, or validate full TOADS inhomogeneity correction.
