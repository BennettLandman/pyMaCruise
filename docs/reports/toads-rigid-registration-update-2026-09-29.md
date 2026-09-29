# P18bd: Java-composed TOADS rigid-registration update

## Scope

The comparison composes the selected VandyCRUISE membership registration
coefficient loop, coefficient normalization, and `registerGradientDescent`
line search. The Java harness uses the supplied `ParametricTransform`,
`RotationMatrix`, and `ImageFunctions` classes for coordinate mapping,
derivatives, and linear interpolation. The synthetic cases are constructed
from deterministic 4×4×4 membership arrays and 5×5×5 atlas priors; no atlas
image collection or subject images are used.

## Evidence

- The Java harness records initial energy, raw and normalized coefficient
  arrays, every trial transform/energy/lambda/improvement decision, final
  transform and energy, lambda, relative change, and rejection count.
- Python matches the Java reference within float32 `rtol=2e-6`, `atol=2e-7`.
- An improving case accepts its first trial; a backtracked case rejects all
  eleven trials and keeps its initial transform.
- The test verifies the selected source's source-faithful coefficient
  normalization call order, including the historical gradient/Hessian
  behavior.
- Focused test: `.venv/bin/python -m pytest -q
  tests/test_toads_rigid_registration_update.py` — **1 passed**.
- Full test: `.venv/bin/python -m pytest -q` — **365 passed in 85.11 seconds**
  on implementation commit
  `82f8bb2d79a08ae3e3846a8670c6af4866269503`. Full output summary is recorded
  in [the full-suite report](full-pytest-2026-09-29.md).

## Provenance and limits

Harness source hashes for the selected atlas method and copied helpers are
embedded in
[`tests/data/synthetic/toads_rigid_registration_update/java_reference.json`](../../tests/data/synthetic/toads_rigid_registration_update/java_reference.json).
The generator is
[`scripts/generate_java_toads_rigid_registration_update_reference.sh`](../../scripts/generate_java_toads_rigid_registration_update_reference.sh).
The harness and fixture are test-only.

This validates one rigid optimizer update against the selected source behavior.
It does not validate physical alignment between the bundled atlas and input
world geometry, full atlas-scale execution, pyramid optimization, deformable
TOADS, the dura stage, or tissue segmentation. It makes no claim of full TOADS
or historical workflow equivalence.
