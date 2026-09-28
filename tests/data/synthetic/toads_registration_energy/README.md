# TOADS registration-energy synthetic reference

This fixture is a deterministic 5×4×3×4 synthetic membership/prior field. It
contains no imaging data. Values are generated in float32 using the formulas
in `synthetic_registration_objective_inputs()` in
`src/pymacruise/toads_registration.py`.

`java_reference.json` is produced by the test-only source-formula excerpt
`tests/java_reference/ToadsRegistrationEnergyReferenceHarness.java`. It
compares the source x/y/z/class accumulation order, ignored background class,
positive-prior filter, class-registration mask, subsample offset, float32
arithmetic, and `cost / norm`. It does not execute the full private atlas
registration method or validate interpolation, derivatives, rigid transforms,
or optimizer convergence.

Regenerate the reference with:

```sh
mkdir -p /tmp/pymacruise-toads-registration-java
javac -d /tmp/pymacruise-toads-registration-java tests/java_reference/ToadsRegistrationEnergyReferenceHarness.java
java -cp /tmp/pymacruise-toads-registration-java ToadsRegistrationEnergyReferenceHarness > tests/data/synthetic/toads_registration_energy/java_reference.json
```

The selected upstream source SHA-256 is
`e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`.
The committed reference JSON SHA-256 is
`e62d1605439f623c3b9ea09ff78b252540eacd341168b508ad8645303cd7d5d6`.
The harness SHA-256 is
`3c2e8bc403a1e04f201cbc29f6aa468175b202c4d2013bc1aefa4429d7bf8b41`.
