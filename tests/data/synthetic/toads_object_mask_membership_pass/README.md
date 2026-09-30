# TOADS OBJECT/MASK membership-pass reference

This synthetic `4×4×4` two-modality, three-class fixture selects OBJECT, MASK,
OBJECT classes in one `computeExactMemberships` pass. It exercises the MASK
binary mismatch energy alongside field/PSF-corrected OBJECT residuals,
source-ordered in-place six-face smoothing, atlas priors, relation factors,
reciprocal normalization, and objective accumulation. The Java harness
checks source SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
All data are formula-generated and synthetic; no image files are used.
Regenerate with `scripts/generate_java_toads_object_mask_membership_pass_reference.sh`.
