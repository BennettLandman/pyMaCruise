# LesionTOADS adaptive field update fixture

`java_reference.json` records source-excerpt outputs for the LesionTOADS
post-thinning degree gate and channel call order. The pinned excerpt is
`LesionToadsInhomogeneityTransitionReferenceHarness.java`.

The pytest test constructs deterministic 7×7×7 float32 attached correction
fields for two channels and mutates their interior through a test refit
callback. This isolates the wrapper scheduling and attached-field identity;
it does not represent measured MRI data or validate the numerical polynomial
fit itself. The source correction class is `InhomogeneityCorrection2`; its
separate exact-fit kernel remains a dependency for a full 13-class integration.
