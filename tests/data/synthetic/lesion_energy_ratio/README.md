# LesionTOADS energy-ratio fixture

This scalar fixture connects the selected initial-thinning objective from
`lesion_initial_thinning_composition/java_reference.json` (P39ak) to the
post-thinning GMflag=true membership objective from
`lesion_post_thinning_gm_membership_update/java_reference.json` (P39am).
The main chain uses prior energy bits `426d9692` and current energy bits
`7fc00000`; the selected float32 wrapper expression produces ratio bits
`ffc00000`.

The synthetic P39am profile has a NaN objective because absent classes have
zero variance. That is the actual selected-Java input to this ratio seam.
Additional scalar controls pin a finite ratio with either operand ordering,
equal energies, `0/0`, and nonzero-over-zero behavior. The harness emits raw
bits for both inputs and each arithmetic intermediate. It also captures the
selected call order:

1. `Eprev = energy`
2. `energy = segmentation.computeMemberships()`
3. `Eratio = Numerics.abs(2.0f * (Eprev - energy) / (Eprev + energy))`

The local Java excerpt of `Numerics.abs(float)` is pinned to the selected
VandyCRUISE `Numerics.java` SHA-256 in the JSON reference and harness. Its
`a > 0 ? a : -a` behavior turns a positive quiet NaN into a negative quiet
NaN and positive zero into negative zero. This matches the Java source and is
different from `Math.abs`/NumPy absolute-value behavior for those bit signs.

This fixture covers only the float32 ratio state update. It does not include
registration, convergence/iteration decisions, growth, classification, or a
complete LesionTOADS iteration. No image is needed because P39ak and P39am
already validate the upstream image-derived objective values independently.
