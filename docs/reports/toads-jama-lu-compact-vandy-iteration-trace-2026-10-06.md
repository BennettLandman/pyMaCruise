# TOADS JAMA LU and compact Vandy iteration trace — 2026-10-06

## Change

`InhomogeneityCorrection2` solves the selected polynomial normal equations through `Jama.Matrix.solve`; square systems use JAMA's `LUDecomposition.solve`. PyMaCRUISE previously delegated these systems to NumPy/LAPACK. Their coefficients could differ by a few float64 ULPs, which changed later float32 field values and eventually caused a compact TOADS loop divergence. The field-fit path now follows JAMA's left-looking LU update order, strict-greater-than partial-pivot selection, and forward/back substitutions.

The standalone degree 1–3 Java fit fixtures continue to use their documented coefficient/field tolerances. This change does not claim universal bitwise coefficient equality across those separate fixtures. On the selected compact Vandy loop trace, however, the resulting corrected fields, memberships, and rigid transforms now match the instrumented Java run exactly.

## Synthetic comparison

- Evidence type: selected VandyCRUISE method/source comparison on one compact synthetic fixture; not whole-workflow equivalence or real-MRI validation.
- Input: committed T1 and matching mask with native shape `(36, 45, 34)`; the selected Java algorithm pads to `(53, 60, 51)` and uses 11 atlas classes.
- Loop settings: one registration level, one registration update per cycle, relation limit 3, and a 12-iteration cap.
- Java references: `MultipleToadDeformableAtlas.java` SHA-256 `e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`; `AlgorithmMultipleToads.java` SHA-256 `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`.
- JAMA reference: `tests/java_reference/upstream/Jama/LUDecomposition.java`, SHA-256 `9b41fefaef78fdc95a2c8fdfebb1aaee8375dcaf6c7b8eddb30acdf7f5a077a6`.

## Results

- All 12 Python membership arrays match the selected Java arrays bit-for-bit (`53 × 60 × 51 × 11` float32 values per iteration).
- All six rigid-transform parameters match their Java float32 bit patterns after each of the 12 cycles.
- The instrumented Java and Python correction fields match exactly for captured iterations 4–12.
- Final native-grid classification matches exactly: zero changed voxels, Dice 1.0 for all 11 labels, and SHA-256 `aff69a659d72b3733edf3d6913c8bd7212d445c977dcc4c2dfe2b2cdb7a2e1c3`.
- Both runs reach the 12-iteration limit without reaching the energy threshold. The Python terminal energy ratio is `0.023335229605436325`; the run must be described as capped, not converged.

The compact case improves the prior iteration-4 divergence to exact equality across the recorded loop. It does not exercise the complete production atlas/workflow configuration, establish natural convergence, or replace the open whole-workflow and real-MRI comparisons.

## Validation

Focused command, run on implementation commit `52b68f5d1b74d723c2a93e76031b16ac2e25e193`:

```bash
./.venv/bin/pytest -q \
  tests/test_toads_inhomogeneity_fit_degree1.py \
  tests/test_toads_inhomogeneity_fit_degree23.py \
  tests/test_toads_inhomogeneity_transition.py \
  tests/test_toads_multichannel_exact_field_refit.py \
  tests/test_toads_vandy_iteration_boundary.py \
  tests/test_base_workflow_restartable.py
```

Result: **16 passed, 0 failed, 0 skipped in 67.51 seconds**. This includes selected degree 1–3 field-fit references, singular-fit rejection, the 12-cycle Java membership/pose/final-label comparison, and restartable synthetic base-workflow checks. No full suite was run, so the full-suite pass count and selected TOADS source-line estimate in `plan.md` remain unchanged.
