# Default TOADS first thinning pass

This analytic 36³, two-channel, 11-class fixture targets the selected
`MedicAlgorithmToads08` default constructor route. The image and initial label
map are stored as RAS `.nii.gz` files with the same anisotropic millimetre
affine. The discrete label field combines linear and cross-axis modular terms
to create interleaved folded interfaces instead of concentric objects.

The Java/Python composition evaluates the selected exact reciprocal-distance
initializer on the source's ten-voxel interior crop, using the fixed 2×11
centroids recorded in `spec.json`. Those full-class memberships feed directly
into the already source-matched first `initThinningLabels` plus
`propagateThinning` pass. The default wrapper's volume ratio is 0.5 and its
thinning priority limit is 3. The selected 18/6 critical-point LUT and heap
ordering are the same assets used by P18az. Relation flags and transformed
shape values are fixed to zero so this fixture isolates the exact-membership
to thinning handoff; registration, relation initialization, field correction,
statistics refresh, next membership update, growing, and output generation are
outside scope.

Regenerate the NIfTI phantom with
`scripts/generate_toads_default_first_thinning_fixture.py`; regenerate the
pinned selected-method JSON with
`scripts/generate_java_toads_default_first_thinning_reference.sh`. The harness
executes source-order exact initialization and thinning formulas, including
the source minimum `max(100, 0.5 × initial class volume)`. It records hashes of
every float32 membership and of the complete queue trace, final segmentation,
skeleton, availability, and residual queue, alongside class-volume bit
patterns. The harness is a selected-method composition reference, not a Java
workflow run or a claim of full TOADS equivalence.

`fixed_inhomogeneity_field.nii.gz` supplies a deterministic non-unit field for
the following fixed-input centroid-method comparison. It uses the selected
T1-MPRAGE descriptor factors (CSF 0.8, GM 1.0, WM 1.0), routes MASK through
its unpowered generic centroid update and the ten optimized classes through
the `PowerTable` kernel, then checks all 11 damped centroids against a Java
source-method harness (`scripts/generate_java_toads_default_exact_centroid_reference.sh`).
The first thinning fixture leaves categorical labels
unchanged, so those labels are the observed post-thinning state. This does not
implement field fitting or later membership updates. The next selected default
stage, `computeExactVariances`, is covered by a separate Java source-method
harness (`scripts/generate_java_toads_default_exact_variance_reference.sh`). It
uses the post-thinning assigned labels, initialized memberships, fixed field,
and Java-verified refreshed centroids. The source scans the full grid for each
assigned class and computes the membership-squared weighted corrected-intensity
standard deviation, including MASK and the ten optimized classes in the default
descriptor. The 11 output float32 bit patterns are checked independently. This
fixed-state composition does not implement subsequent exact-membership updates,
registration, stop checks, growing, or final label classification.

`membership_transformed_shapes.nii.gz` and `membership_closest_labels.nii.gz`
pin a folded shape/relation state for the next exact default membership slice;
their formulas and all stage statistics/scalars are recorded in
`membership_stage_spec.json`. The composed adapter chains the two-channel
initial exact membership/thinning state, exact centroid refresh, assigned-class
standard deviations, centroid similarity (with explicit image-variance
weights), then the existing 11-class exact membership sweep. The independent
Java composition harness checks centroid and deviation bits, similarity/data
energy/denominator/membership SHA-256 values over every float bit, and the
objective bits. This fixed-input selected-method comparison leaves the
wrapper's image-wide variance estimation, field refit, atlas registration,
later iterations/convergence, growing, and final output classification out of
scope. The Java composition harness relies on the separately pinned first
thinning reference for this fixture's unchanged categorical labels.
