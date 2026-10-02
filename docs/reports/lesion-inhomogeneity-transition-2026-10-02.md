# LesionTOADS adaptive inhomogeneity transition

This slice ports the selected post-thinning control branch in
`AlgorithmLesionToads.runAlgorithm`: after `propagateThinning`, the wrapper
compares the previous float32 energy ratio with `lowEratio`, where the
constructor sets `lowEratio = (float)Math.sqrt(maxDistance)`. It increments
the current polynomial degree by exactly one when the ratio is strictly below
that threshold, then calls `computeCorrectionField(deg)` for each image
channel in order. Degree zero is still dispatched; the selected Java
`InhomogeneityCorrection2.computeExact3DPolynomialField(0)` immediately
returns.

`src/pymacruise/lesion_inhomogeneity_transition.py` preserves the float32
threshold conversion, strict comparison, cap, ordered per-channel calls, and
the identity of each mutable attached field. The polynomial fit is injected
as a channel callback, because the current exact Python correction fitter is
limited to the 11-class base-TOADS descriptor while LesionTOADS uses a
13-class profile. This report therefore records a **branch-level source
match**, not an end-to-end field-fit or LesionTOADS workflow equivalence.

## Evidence

- Focused command: `./.venv/bin/pytest -q tests/test_lesion_inhomogeneity_transition.py`
- Result: **2 passed in 0.91 seconds**.
- The test compiles and runs
  `tests/java_reference/LesionToadsInhomogeneityTransitionReferenceHarness.java`
  and compares its output with
  `tests/data/synthetic/lesion_inhomogeneity_transition/java_reference.json`.
- Harness-pinned source hashes: `AlgorithmLesionToads.java`
  `556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`;
  `InhomogeneityCorrection2.java`
  `6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8`.
- Harness SHA256:
  `9a4e519a5d98c397c141ed36bcd35788d5527a1ba46c30c65b469dff7f85c84c`.

Cases cover float32 values immediately below, at, and above the threshold, a
degree increment, the configured degree cap, degree-zero dispatch, and a NaN
energy ratio (which does not satisfy Java's `<` comparison). The test refit
callback only verifies channel order and field mutation; numerical polynomial
refit results remain unvalidated for the 13-class LesionTOADS profile.
