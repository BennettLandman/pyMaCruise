# Optimized TOADS image-variance phantom

This source-formula reference reuses the 7×6×5×3 RAS, 1 mm isotropic
synthetic image and multiplicative-field NIfTI files in
`../toads_optimized_likelihood/`, plus the first three membership channels
from `../toads_direct_membership/`. No duplicate image volumes are stored.
The harness defines three optimized tissue classes and uses the CSF/GM/WM
centroids from the selected TOADS initialization fixture.

`java_reference.json` records modality-normalized inverse variance weights
and the unnormalized inverse-variance square roots with the field off/on. It
is generated from the selected `computeExactImageVariance` optimized-class
formula excerpt. The selected `MultipleToadSegmentation.java` SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
`INF=1e10f` behavior is preserved for a modality with zero numerator; in
multi-modality zero-residual cases those values are still normalized by their
finite sum, matching the source.

Regenerate the reference with
`scripts/generate_java_toads_optimized_variance_reference.sh`. These are
synthetic source-formula checks, not a full `MultipleToadSegmentation` run.
