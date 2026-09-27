# Surface extraction phantoms

These committed `.nii.gz` level sets are deterministic synthetic inputs for
P30. `sphere.nii.gz` has analytic RAS+ radius 5.5 mm. `ellipsoid.nii.gz` has
analytic RAS+ radii 5.5, 4.5, and 6.0 mm. Both use a 31³ grid at 0.5 mm
spacing and exact signed implicit fields; they contain no subject data.

Regenerate with `python tests/data/synthetic/surfaces/make_phantoms.py` from
the repository root. Verify artifacts using `shasum -a 256 -c SHA256SUMS`.
