# LesionTOADS unfielded data-energy formula reference

The committed Java fixture has nine deterministic three-channel cases for
the selected unfielded `OBJECT`, `LESION`, and `BLACKHOLE` branches in
`LesionToadSegmentation.computeMemberships`. It records intensity and
centroid inputs, nonuniform channel weights, and the raw float32 class,
lesion, and black-hole energy accumulators as IEEE bit patterns.

Regenerate the reference with
`scripts/generate_java_lesion_toads_data_terms_reference.sh`. These cases
validate only unfielded squared-residual data terms before inversion,
PowerTable weighting, spatial smoothing, priors, or membership normalization.
