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
The separate coefficient call-order excerpt is
`tests/java_reference/ToadsRegistrationCoefficientReferenceHarness.java`,
SHA-256 `cc62e5e25bd3297a0d03a981c951cea0f7bb137810fe0791b2952114c39af1a5`.
Its output fixture is `java_coefficient_reference.json`, SHA-256
`a0288b795161d02d42e83cda857ce61355924fbb597e172d7c0b8781bafba368`.
The prior-derivative excerpt is
`tests/java_reference/ToadsPriorDerivativeReferenceHarness.java`, SHA-256
`127099ca604f07537992594a9dbe871ce0efe0fdf3a3106277ae91c72f85a5b0`; its
four synthetic reference points are in
`java_prior_derivative_reference.json` (SHA-256
`957162ea13b990895098c61a81876a87884c26d22834e570f51476d14718d849`).
The rigid-coordinate Jacobian excerpt is
`tests/java_reference/ToadsRigidJacobianReferenceHarness.java`, SHA-256
`e04a9319eae051362547f0461100512a0bcf58e0d7f70763c300238f939bbcc4`; its
synthetic reference is `java_rigid_jacobian_reference.json` (SHA-256
`db9450821c6142e7383b94b27a158224efc6531f4980a04bd89a0c0251d542a3`).
The one-voxel coefficient excerpt is
`tests/java_reference/ToadsRegistrationCoefficientAccumulationHarness.java`,
SHA-256 `20ac7238c1eaac7c6a9b72a2a974cca58582c0eaf2ab8fd7dead86755d7aa2d9`;
its output is `java_coefficient_accumulation_reference.json` (SHA-256
`fbef8094da69fe49d512afd143810dddb18e6400901a57629da13e7b654afa60`).
