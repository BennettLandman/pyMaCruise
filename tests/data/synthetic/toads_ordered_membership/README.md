# TOADS ordered membership-sweep reference

This comparison reuses the deterministic committed RAS `.nii.gz` inputs from
`toads_direct_membership/` and `toads_shape_relations/`. The Java formula
harness outputs a 7×6×5×4 updated membership volume and one cumulative
float32 objective scalar for a single small-grid sweep.

The harness visits x, then y, then z, then class, using the in-place
membership array for six-face smoothing. It combines the existing synthetic
data energies, smoothing coefficient/priors, shape relation factors, source
`ZERO`/`INF` rule, class normalization, and source objective accumulation.
Source file hashes and the exact selected parameter values are included in
`java_reference.json`.

Regenerate the reference with
`scripts/generate_java_toads_ordered_membership_reference.sh`. This is a
source-formula excerpt for one sweep, with precomputed data terms and no
topology admissibility check; it is not full TOADS segmentation equivalence.
