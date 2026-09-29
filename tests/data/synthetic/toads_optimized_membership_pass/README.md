# Optimized-tissue TOADS membership-pass phantom

This composition reference reuses the committed 7×6×5×3 RAS 1 mm image and
multiplicative field volumes in `../toads_optimized_likelihood/`, the first
three membership channels in `../toads_direct_membership/`, the categorical
assignment map in `../toads_optimized_class_variance/`, and the first three
shape and closest-label channels in `../toads_shape_relations/`. No duplicate
image volumes are stored.

The Java harness composes source-formula global image variance, assigned
class deviations, centroid similarities, PSF ridge bias, actual
`PowerTable`-based optimized likelihoods, atlas inverse-shape-prior energy,
six-face smoothing, relation weighting, reciprocal conversion, and the
ordered normalized update. It generates separate JSON references with the
field enabled and disabled. The selected
`MultipleToadSegmentation.java` SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
Regenerate with
`scripts/generate_java_toads_optimized_membership_pass_reference.sh`.

This validates one small optimized-tissue membership pass only. It is not a
complete TOADS segmentation or topology-constrained evolution run.
