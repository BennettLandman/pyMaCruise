# TOADS JAMA LU and compact Vandy convergence trace — 2026-10-06

## Change

`InhomogeneityCorrection2` solves the selected polynomial normal equations through `Jama.Matrix.solve`; square systems use JAMA's `LUDecomposition.solve`. PyMaCRUISE previously delegated these systems to NumPy/LAPACK. Their coefficients could differ by a few float64 ULPs, which changed later float32 field values and eventually caused a compact TOADS loop divergence. The field-fit path now follows JAMA's left-looking LU update order, strict-greater-than partial-pivot selection, and forward/back substitutions.

The standalone degree 1–3 Java fit fixtures continue to use their documented coefficient/field tolerances. This change does not claim universal bitwise coefficient equality across those separate fixtures. On the selected compact Vandy loop trace, however, the resulting corrected fields, memberships, and rigid transforms now match the instrumented Java run exactly.

## Synthetic comparison

- Evidence type: selected VandyCRUISE method/source comparison on one compact synthetic fixture; not whole-workflow equivalence or real-MRI validation.
- Input: committed T1 and matching mask with native shape `(36, 45, 34)`; the selected Java algorithm pads to `(53, 60, 51)` and uses 11 atlas classes.
- Loop settings: one registration level, one registration update per cycle, relation limit 3, and a 50-iteration maximum. Both implementations meet the source threshold before the cap.
- Java references: `MultipleToadDeformableAtlas.java` SHA-256 `e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`; `AlgorithmMultipleToads.java` SHA-256 `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`.
- JAMA reference: `tests/java_reference/upstream/Jama/LUDecomposition.java`, SHA-256 `9b41fefaef78fdc95a2c8fdfebb1aaee8375dcaf6c7b8eddb30acdf7f5a077a6`.

## Results

- All 19 Python membership states match the selected Java C-order float32 hashes exactly (`53 × 60 × 51 × 11` values per iteration); complete Java membership arrays are also committed for iterations 1–12.
- All six rigid-transform parameters match their Java float32 bit patterns after each of the 19 cycles.
- The Java and Python final terminal ratio both round to float32 `9.45949e-5`, below the configured `1e-4` threshold; both stop naturally at iteration 19, before the 50-iteration cap.
- Final native-grid classification matches exactly: zero changed voxels, Dice 1.0 for all 11 labels, and SHA-256 `170709d51a984c74a6f1a400c7ca156045e62f54a94f9d7131adef628206bfb6`.

The JAMA solver change improves the prior iteration-4 divergence to an exact membership and transform match through natural convergence on this compact selected-source case. It does not establish full MaCRUISE workflow equivalence or real-MRI acceptance.

## Validation

Focused command, run on validation commit `625e86b5be8cf883a3b475f2c2ba9421772b15c1` (solver implementation commit `52b68f5d1b74d723c2a93e76031b16ac2e25e193`):

```bash
./.venv/bin/pytest -q \
  tests/test_toads_inhomogeneity_fit_degree1.py \
  tests/test_toads_inhomogeneity_fit_degree23.py \
  tests/test_toads_inhomogeneity_transition.py \
  tests/test_toads_multichannel_exact_field_refit.py \
  tests/test_toads_vandy_iteration_boundary.py \
  tests/test_base_workflow_restartable.py
```

Result: **16 passed, 0 failed, 0 skipped in 114.13 seconds**. This includes selected degree 1–3 field-fit references, singular-fit rejection, the 19-cycle Java membership/pose/final-label comparison through natural convergence, and restartable synthetic base-workflow checks. The latest full suite on tree `2b1d0c51888725f3d6a79138ac777bab06bb0b0d` passed **1,218 tests, 0 failed, 0 skipped, 11 expected RAS warnings in 1,452.34 seconds**, with **10.179 GB** peak sampled RSS; see [full-suite evidence](full-pytest-2026-10-06-jama-lu.md). The suite result covers the solver implementation; the 19-cycle source-reference extension itself received the focused validation above. It adds zero Vandy source lines, so the source-line estimate remains unchanged.
