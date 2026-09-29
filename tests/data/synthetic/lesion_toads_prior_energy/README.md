# LesionTOADS atlas-prior energy reference

Five selected-Java one-voxel cases exercise the `priorCoefficient` class,
lesion, and black-hole shape/similarity accumulators. They cover ordinary,
lesion, optimized-lesion, and black-hole target classes; the `m != k`
exclusion; optional lesion-class terms; the extra black-hole similarity term
when the selected lesion class is a black hole; zero priors/shapes; and
float32-sensitive values.

The fixture starts from provided shape and similarity vectors. It does not
construct similarity maps, perform registration, or compute an entire image
volume. Regenerate it with
`scripts/generate_java_lesion_toads_prior_energy_reference.sh`.
