# TOADS degree-one approximate field-fit fixture

The 9 x 9 x 9 RAS phantom uses an anisotropic affine with voxel sizes 0.85,
1.15, and 1.55 mm. It contains a positive scalar image generated from a known
linear field, eleven centroids, three membership slots per voxel, and a
nonempty segmentation mask. Candidate slot zero is always atlas class 0, the
background class that the selected constructor identifies from its topology;
the fit skips it. Positive classes have a shared centroid so the fitted degree
one model is identifiable and can be checked against a known field.

The pinned formula harness covers fixed-degree
`computeApprox3DPolynomialField(1)`. Although the constructor requests IMAGE,
the selected default atlas contains `optimized_*` classes and changes the
effective factor to OPTIMIZED_IMAGE. With null optional category exponents,
the three optimized category PowerTables use p=1. The fit also uses
GLOBAL/CHEBYSHEV/FIXED/EQUAL_VARIANCE, q=2 settings. It performs the source
stride-3 fit, full-volume system scaling, a square Jama-style LU solve,
full-grid evaluation, and correction-mask mean centering. Source hashes for
`InhomogeneityCorrection2`, `JamaMatrix`, `Jama.Matrix`, `LUDecomposition`, and
`PowerTable` are embedded in `java_reference.json`; the reference harness
SHA-256 is
`80b40a15ee803fba404aba312307b15a6b7c7d93c230d2c6432d7eab770c6910`.

The normal matrix and right-hand side match the Java float/double accumulation
exactly. Coefficients differ by at most `2.92e-16` because NumPy and the
Jama-style reference LU solve have slightly different double-precision
rounding. The mean-normalized float32 field matches the Java output exactly.
Its float64 measured mean is within `6.56e-7` of one after the source's
sequential float32 mean and offset operations.

This validates the fixed degree-one fit only. Degree two/three, singular
behavior beyond explicit rejection, the degree ramp, segmentation field
insertion, later iterations, spline mode, and full TOADS equivalence remain
outside scope.
