# TOADS OBJECT/MASK/OUTLIER membership pass

This deterministic 4 x 4 x 4 synthetic case combines the selected OBJECT
residual, MASK mismatch, and paired OUTLIER energy terms with field correction,
PSF bias, six-face smoothing, atlas inverse-shape prior, relation factors,
class normalization, outlier-map normalization, and objective accumulation.
It reproduces the selected Java dangling-else behavior: smoothing is added to
the ordinary OUTLIER energy but skipped for the auxiliary energy; the separate
prior term updates both. A fixed interior probe in the Java reference records
the raw energies, smoothing increment, prior increment, and resulting
auxiliary energy as float bits.
It has one modality and three classes and contains no subject data.

Regenerate with
`scripts/generate_java_toads_object_mask_outlier_membership_pass_reference.sh`.
The generator checks the selected
`MultipleToadSegmentation.java` SHA-256 before compiling the Java formula
harness. This compares one selected synthetic membership pass; it does not
validate TOADS topology, deformation, or outer-loop convergence.
