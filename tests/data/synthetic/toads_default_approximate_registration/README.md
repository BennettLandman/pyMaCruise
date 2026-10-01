# P18czp/P18czq: default TOADS no-atlas initialization and initial pyramid

This asymmetric synthetic case validates the selected default 11-class
approximate initialization and a method-level initial rigid registration
composition, with source call order taken from AlgorithmMultipleToads.runAlgorithm():

1. no-atlas MultipleToadPreprocess.initialMemberships(...)
2. atlas.setImageInfo(...)
3. atlas.initShapeRegistration(...)
4. set the atlas transform from preprocessing
5. atlas.updateShapeRegistration(...)
6. atlas.registerShapesPyramid(), dispatched to registerShapesApproxPyramid()

P18czq composes the P18czp intensity-only initializer with the selected rigid
atlas-registration pyramid. The Java harness pins the wrapper's source call
order and runs the selected membership-pooling, atlas-shape downsampling,
registration coefficient, and line-search excerpts. It starts from a known,
nonzero synthetic pre-pyramid rigid transform, modeling the value handed from
`preprocess.getTransform()` to the atlas. It does not recompute that transform
from source image orientation or crop state, nor execute the Java wrapper's
`initShapeRegistration`/`updateShapeRegistration` setup methods. The Python
composition test calls `initialize_and_register_toads_default_approximate_atlas`
and checks that its initialization arrays and per-scale transform handoff match
the selected reference. This path does not call ANTsPy; the internal TOADS
objective and optimizer use the selected Python kernels. ANTsPy is used for
other workflow registration boundaries, which have separate transform and
interpolation validation. Neither this fixture nor the selected-method match
establishes full TOADS or end-to-end registration equivalence.

P18czp covers only the no-atlas approximate overload selected by
AlgorithmMultipleToads for its initial pyramid memberships. In particular,
the exercised method is `MultipleToadPreprocess.initialApproxMemberships`
through `initialMemberships(memberships,best,approx,centroids,classes)`; this
intensity-only initialization does not inspect atlas shape maps. The dedicated
AlgorithmLesionToads wrapper instead invokes
`LesionToadPreprocess.initialMemberships(centroids,classes)`, a separate
full-class method with no `approx` argument. That LesionTOADS initializer is
not covered by P18czp; its implementation and evidence require a separate plan
row/slice.

The separate atlas-aware initial-membership overload used for later
segmentation memberships is outside this fixture's call path. This validates
the selected initialization/pyramid harness path only; it does not establish
full TOADS workflow equivalence or segmentation accuracy.

The P18czp Java reference also pins the zero-residual tie, below-`1e-20`
residual cutoff, and repeated-centroid multiplicity cases. It confirms that
the overload applies neither an atlas shape penalty nor a forced candidate
slot. In the asymmetric volume comparison, all membership bits and candidate
indices are compared, including the first and last boundary voxels.

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
  It also pins the compressed subject and atlas-shape inputs by SHA-256.
  The test also pins the harness SHA-256:
  9e52b8511b2fe1f874d15d60691873ee064169342d080478d446c3911b894827.

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
