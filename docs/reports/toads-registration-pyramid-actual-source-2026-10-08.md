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
source method's final six float32 transform parameters and final line-search
lambda exactly.

The coarse level has zero gradient and Hessian for this small phantom. The
source takes its 11 rejected NaN trials before advancing to scale 2. Python
preserves the unchanged coarse transform and rejection count, but currently
reports a NaN coarse-level energy where the Java path produces a finite
decreased energy while evaluating NaN coordinates. Thus this evidence does
not yet validate per-level energies, relative-change values, or complete
optimizer state. P18ck remains Started until those intermediate states are
reconciled and the broader suite is rerun.

## Validation

- Direct production-method comparison:
  `pytest -q tests/test_toads_registration_sampling.py::test_pyramid_registration_matches_actual_vandycruise_entry_point` —
  **1 passed**.
- Focused registration module:
  `pytest -q tests/test_toads_registration_sampling.py` — **22 passed**.
- Alignment-startup regressions:
  `pytest -q tests/test_toads_alignment_startup.py` — **6 passed**.
- These focused checks are synthetic source comparisons. They do not validate
  full TOADS segmentation, the full MaCRUISE workflow, or real MRI data.
- The most recent full suite predates this change; no full-suite claim is made
  for this report.
