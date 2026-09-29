# Ordered TOADS prior-energy sweep phantom

This composition test reuses the committed 7×6×5×4 RAS 1 mm NIfTI fixtures
for synthetic data energy and initial memberships from
`../toads_direct_membership/`, and transformed shape values/closest labels
from `../toads_shape_relations/`. It adds four atlas shape-prior weights and
a symmetric 4×4 centroid-similarity matrix, computed from the selected
source formula with fixed synthetic centroids, standard deviations, modality
weights, and `priorScale=0.05f`.

`java_reference.json` records the final membership volume and cumulative
objective after one ordered in-place sweep with six-face smoothing and atlas
inverse-shape-prior energy. The harness composes the selected source formulas
for data energy, smoothing, shape prior, relation factors, reciprocal
conversion, normalization, and objective accumulation. The selected
`MultipleToadSegmentation.java` SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
Regenerate it with
`scripts/generate_java_toads_ordered_prior_sweep_reference.sh`.

This is a compact formula-composition case, not a full TOADS segmentation or
topology-constrained evolution.
