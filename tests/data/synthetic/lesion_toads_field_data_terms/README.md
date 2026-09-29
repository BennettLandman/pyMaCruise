# LesionTOADS field-corrected data-energy formula reference

The committed Java fixture contains nine deterministic three-channel cases
for the selected `useField` `OBJECT`, `LESION`, and `BLACKHOLE` data-term
branches. Nonuniform field values multiply each image channel before class,
lesion, and black-hole centroid residuals. Inputs and raw float32 accumulators
are recorded as bit patterns.

Regenerate with `scripts/generate_java_lesion_toads_field_data_terms_reference.sh`.
This validates field multiplication and energy accumulation only; special
PowerTable handling, smoothing, atlas priors, inversion, and membership
normalization are excluded.
