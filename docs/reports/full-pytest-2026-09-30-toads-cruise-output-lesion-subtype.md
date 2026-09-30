# Full-suite evidence: TOADS CRUISE output and LesionTOADS subtype maps

Date: 2026-09-30

## Full suite

- Command: `.venv/bin/python -m pytest -q`
- Result: **701 passed, 2 expected RAS reorientation warnings, 168.47 seconds**.
- Tested implementation commit: `a6cc677e3228cd9af7bfa945e55f90707b601118` (`Add TOADS outputs and lesion subtype adapters`).
- Both warnings are intentional: the Dura topology adapter receives LPS-oriented topology and paint-mask inputs and verifies normalization to RAS.

Focused adapters: `.venv/bin/python -m pytest -q tests/test_toads_cruise_outputs.py tests/test_lesion_subtype_map.py` — **8 passed in 0.77 seconds**.

## P18cv: base-TOADS CRUISE output maps

`toads_cruise_outputs.py` matches the selected `MultipleToadSegmentation.generateCruiseOutputs()` method excerpt for CSF, GM, WM-fill, and hard-WM maps in exact and approximate modes. The fixtures use asymmetric RAS `.nii.gz` volumes and a pinned Java reference harness. Tests cover aliases, duplicate/absent labels, candidate-slot ordering, exclusions, float32 values, and hard-WM behavior. Wrapper uncrop/fill remains the separately validated P18cn operation. This validates the output adapter; it does not establish TOADS evolution or full CRUISE equivalence.

## P39zzzzi: stopping-step lesion/black-hole subtype map

`lesion_subtype_map.py` matches the selected `updateLesionClassifcation()` excerpt, including clearing stale output, strict membership thresholds, lesion-favoring subtype ties, lesions-disabled behavior, and qualifying-voxel count. Synthetic cases are true-RAS `.nii.gz` and compared against the pinned Java harness. This validates one stopping-step map, not full LesionTOADS segmentation.

## Commit provenance

The implementation commit above is the code state exercised by pytest. The plan-only commit adding the next P18cw task occurred while the run was in progress and did not change Python or test files.
