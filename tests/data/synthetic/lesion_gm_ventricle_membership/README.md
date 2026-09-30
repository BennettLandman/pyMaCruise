# LesionTOADS GM ventricle membership excerpt

This scalar synthetic fixture isolates the ordinary-class branch in
`LesionToadSegmentation.computeMemberships` that attenuates the exact class
names `Cerebrum-GM` and `GM` by `VentDist / MaxVentDist` when a lesion exists,
the maximum is nonzero, and the signed-byte distance is not `-1`. It also
includes the source's second-pass positive-denominator normalization.

The Java harness compares DIRECT and PRIOR paths, both exact names and
near-miss names, zero distance, invalid distance, zero maximum, no lesion,
float32-sensitive arithmetic order, and a nonpositive denominator. Input
`num` terms stand in for values calculated by earlier `computeMemberships`
branches. The fixture does not implement those upstream terms, auxiliary
lesion/black-hole object types, or the complete LesionTOADS segmentation.

Regenerate the reference from the repository root:

```sh
mkdir -p /tmp/pymacruise-gm-ventricle-reference
javac -d /tmp/pymacruise-gm-ventricle-reference tests/java_reference/LesionGmVentricleMembershipReferenceHarness.java
java -cp /tmp/pymacruise-gm-ventricle-reference LesionGmVentricleMembershipReferenceHarness \
  > tests/data/synthetic/lesion_gm_ventricle_membership/java_reference.json
```

The fixture records Java `Float.floatToIntBits` values, so tests check float32
operation order exactly rather than applying a tolerance.
