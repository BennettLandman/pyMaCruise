# TOADS six-face smoothing phantom

This deterministic fixture contains a 7×6×5×4 float32 membership volume in a
1 mm isotropic RAS `.nii.gz` grid. Voxel/class values are defined by
`scripts/generate_toads_smoothing_fixture.py`; the asymmetric values make
axis traversal and class exclusion observable.

`java_reference.json` stores the per-candidate-class smoothing increments to
the reciprocal-energy `num` and objective `energy` from the selected
`MultipleToadSegmentation.computeExactMemberships` six-face smoothing formula.
The Java harness keeps the source X+, X−, Y+, Y−, Z+, Z− and class-loop order,
uses a smoothing coefficient of 0.017, and source smoothing priors
`[0.5, 0.8, 1.1, 0.25]`. The source file SHA-256 is recorded in the JSON.
Border values are zero because the source update only iterates interior voxels.

Regenerate the NIfTI with `scripts/generate_toads_smoothing_fixture.py` and the
Java formula reference with `scripts/generate_java_toads_smoothing_reference.sh`.
