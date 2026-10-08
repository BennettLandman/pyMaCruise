# TOADS pyramid registration: actual VandyCRUISE source check

## Finding

The prior P18ck comparison validated a helper-level transcription, not the
production `MultipleToadDeformableAtlas.registerShapesPyramid()` method. The
transcription used every coarse-grid voxel and the gradient convention from a
different `MultipleToadAtlas` method. The production VandyCRUISE path sets
`subsample = 2` in `updateShapeRegistration()` and its exact coefficient
method divides the raw gradient and Hessian by the prior norm once each.

The Python pyramid wrapper now defaults to the source's offset-zero,
step-two scan and uses the gradient convention in
`MultipleToadDeformableAtlas`. The lower-level `source_faithful` switch is
retained for the separate parent-class behavior; the pyramid wrapper chooses
the deformation-class behavior by default.

## Direct source harness

`tests/java_reference/ToadsPyramidActualRegistrationHarness.java` calls the
supplied source's actual `updateShapeRegistration(..., "rigid", ...)` and
`registerShapesPyramid()` entry points. Its launcher checks hashes for the
selected VandyCRUISE Java classes before compiling. On the committed
asymmetric synthetic RAS memberships and atlas shapes, Python matches the
source method's final six float32 transform parameters, final line-search
lambda, and final energy exactly. At scale 4, it also matches the source
trace's initial energy, final rejected NaN-trial energy, relative change,
unchanged transform handoff, and 11 rejected trials. At scale 2, the harness
counts the same trial and accepted updates as the Python iteration history.

The coarse level has zero gradient and Hessian for this small phantom. The
source produces NaN parameter trials, but its comparison-based `Numerics`
clamping maps those coordinates to a finite edge interpolation. Python now
reproduces that behavior for rejected energy evaluations. This is a synthetic
actual-method comparison; it does not establish full TOADS segmentation or
MaCRUISE workflow equivalence.

## Validation

- Direct production-method comparison:
  `pytest -q tests/test_toads_registration_sampling.py::test_pyramid_registration_matches_actual_vandycruise_entry_point` —
  **1 passed**.
- Focused registration module:
- Combined registration and alignment-startup modules:
  `pytest -q tests/test_toads_registration_sampling.py tests/test_toads_alignment_startup.py` —
  **28 passed**.
- No full-suite claim is made yet for the latest source-semantics change.
- These focused checks are synthetic source comparisons. They do not validate
  full TOADS segmentation, the full MaCRUISE workflow, or real MRI data.
- The full suite is pending after this change.
