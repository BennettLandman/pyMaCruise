# P18czp/P18czq: default TOADS no-atlas initialization and initial pyramid

This asymmetric synthetic case validates the selected default 11-class
approximate initialization and initial rigid registration call sequence from
AlgorithmMultipleToads.runAlgorithm():

1. no-atlas MultipleToadPreprocess.initialMemberships(...)
2. atlas.setImageInfo(...)
3. atlas.initShapeRegistration(...)
4. set the atlas transform from preprocessing
5. atlas.updateShapeRegistration(...)
6. atlas.registerShapesPyramid(), dispatched to registerShapesApproxPyramid()

The separate atlas-aware initial-membership overload used for later
segmentation memberships is outside this fixture's call path. This validates
the selected initialization/pyramid harness path only; it does not establish
full TOADS workflow equivalence or segmentation accuracy.

## Phantom and reference

- images.nii.gz: two-channel synthetic image, 9×8×7 voxels.
- atlas_shapes.nii.gz: 11 descriptor-ordered maps on a distinct 11×10×9 grid.
- Both NIfTI inputs use RAS axes and 0.9×1.2×1.45 mm spacing.
- java_reference.json records the no-atlas initialization arrays, top-k
  candidates, per-scale pooled arrays, registration state, and line-search
  trace. The method excerpt is
  tests/java_reference/ToadsDefaultApproximateRegistrationReferenceHarness.java.
  The Java helper classes used by the harness are preserved in
  tests/java_reference/upstream/.
- The reference pins AlgorithmMultipleToads.java, MultipleToadPreprocess.java,
  MultipleToadDeformableAtlas.java, and the selected default atlas descriptor.
  The test also pins the harness SHA-256:
  973592718d4d55e20aed157909e5deb566c5cf4a3a552e76b6eddaf1a11abd0e.

## Measured comparison

Scale 4 is bit-exact for pooled memberships/candidate indices, registration
coefficients, every trial parameter/energy, accepted or rejected decisions,
and optimizer state.

At scale 2, the first coefficient evaluation has the only coefficient
differences: the gradient maximum is 23 float32 ULPs (absolute maximum
2.3283064e-10 at magnitude 1.55656e-5); the Hessian maximum is 1 ULP
(absolute maximum 4.656613e-10). The first lambda=1 trial is rejected. Its
energy differs by 3 float32 ULPs (absolute difference 1.1175871e-8), and its
parameter vector differs by at most 7.6293945e-6. The next trial is accepted.
The accepted/final energy, lambda, relative-change value, accepted flag,
rejected count, and stop reason match the Java reference exactly. The accepted
parameter vector and final transform differ by at most 3.8146973e-6. The
second scale-2 coefficient pass and all trial energies are bit-exact; its
parameters carry forward the bounded accepted-state difference.

Only these documented scale-2 coefficient, rejected-trial-energy, and
transform comparisons use tolerances. The test enforces exact comparisons for
scale 4 and for all scalar optimizer decisions and state values that match.

## Regeneration and focused validation

From pymacruise:

    PYTHONPATH=src .venv/bin/python tests/data/synthetic/toads_default_approximate_registration/generate_fixture.py
    tests/java_reference/generate_toads_default_approximate_registration_reference.sh
    PYTHONPATH=src .venv/bin/pytest -q tests/test_toads_default_approximate_initial_membership.py tests/test_toads_default_approximate_registration.py
