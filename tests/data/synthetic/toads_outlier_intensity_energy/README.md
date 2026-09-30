# TOADS OUTLIER intensity-energy reference

This deterministic `3×2×2` synthetic array fixture compares the selected
`MultipleToadSegmentation.computeExactMemberships` OUTLIER data branch from
VandyCRUISE source SHA-256
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
Four cases exercise raw and field-corrected residuals, optional PSF bias, and
two finite scalar outlier likelihoods. The Java harness records exact float32
bits for both residual and outlier energies. No image files or subject data
are used; arrays are formula-generated inside the committed Java harness.

Regenerate with `scripts/generate_java_toads_outlier_intensity_reference.sh`.
