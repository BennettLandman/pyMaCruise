# Folded-cortex TOADS registration pyramid source comparison

Date: 2026-10-08  
Scope: selected VandyCRUISE rigid pyramid method on a larger synthetic folded cortex

## Fixture and command

The committed `(48,48,48,4)` float32 RAS fixture has soft background, white-matter, folded gray-matter, and CSF channels. Its atlas priors are shifted 1.5 voxels from the memberships, so registration makes accepted updates at all three pyramid scales. Generate it with:

```bash
./.venv/bin/python tests/data/synthetic/toads_registration_pyramid_folded/generate_fixture.py
```

Focused validation:

```bash
./.venv/bin/python -m pytest -q tests/test_toads_registration_pyramid_folded.py
```

The test directly runs `ToadsPyramidActualRegistrationHarness`, which calls the supplied VandyCRUISE `MultipleToadDeformableAtlas.updateShapeRegistration()` and `registerShapesPyramid()` methods, then compares that trace with the Python implementation using the safe default that rejects nonfinite trial transforms.

## Result

- **1 passed in 5.03s**.
- Java and Python final six float32 transform values match bit-for-bit.
- At source levels 4, 3, and 2 (Python block scales 8, 4, and 2), the number of line-search trials and accepted updates matches the Java trace; each scale makes at least one accepted update and all returned transforms are finite.
- The last accepted energy and terminal lambda match the Java reference exactly. The Java object's lambda includes the final rejected trial, while the accepted transform and energy remain unchanged.
- This adds a larger folded case beside the existing small degenerate case, which validates the source-faithful nonfinite-trial switch.

This is selected-source synthetic method evidence only. It does not validate atlas-to-subject physical alignment, full TOADS initialization/iteration, the composite MaCRUISE workflow, or real MRI. Selected VandyCRUISE source-line estimate delta: **0**.
