# TOADS relation-factor fixture

This fixture records selected Java helper outputs for the three scalar factors
used by `MultipleToadSegmentation.computeExactMemberships`: the signed-byte
boundary-distance factor, shape-prior cap, and bounded neighbor-prior factor.
It covers boundary values from -128 to 127, four positive float32 marching
limits, clamp thresholds, negative/zero/positive shape values, and five
float32 prior factors. There is no medical image data.

Regenerate with `scripts/generate_java_toads_relation_factors_reference.sh`.
Source SHA-256 values are recorded in the fixture for
`MultipleToadSegmentation.java` and `Numerics.java`.
