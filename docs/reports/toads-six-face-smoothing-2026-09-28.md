# TOADS six-face smoothing comparison

**Date:** 2026-09-28  
**Roadmap:** P18ag complete  
**Validated implementation/fixture commit:** `d8a02be707cc0a14ed1876e6bb7b3476aba71a43`

## Selected source behavior

`src/pymacruise/toads_math.py::toads_six_neighbor_smoothing_terms` ports the
spatial smoothing contribution in VandyCRUISE
`MultipleToadSegmentation.computeExactMemberships`. For each interior voxel
and candidate class `k`, it visits the six face neighbors in this order:
X+, X−, Y+, Y−, Z+, Z−. At each neighbor it accumulates the squared
float32 memberships of every class except `k`. It then returns the source
increments

- `num += smoothingPrior[k] * smoothing * neighborSum`
- `energy += 0.5f * smoothingPrior[k] * smoothing * neighborSum`

with float32 operation and accumulation order. Source update loops leave all
outer voxels untouched. The source guard only applies smoothing when its
coefficient is positive, so nonpositive values return zero increments.

The selected `MultipleToadSegmentation.java` source SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
`ToadsSmoothingReferenceHarness.java` transcribes the selected six-neighbor
loops and term additions. Its fixture records four class priors and the source
coefficient used in this slice.

## Synthetic fixture and results

The committed `tests/data/synthetic/toads_smoothing/memberships.nii.gz` is a
7×6×5×4 float32 RAS volume with a 1 mm isotropic affine. Asymmetric synthetic
memberships make axis ordering and candidate-class exclusion observable. All
840 values are checked against the deterministic fixture recipe.

Both the `num` and `energy` increments match the Java formula reference
exactly for every voxel and class. Border voxels remain zero. Six separate
one-neighbor perturbations prove that each face contributes; changing a
neighbor membership for class 1 leaves candidate class 1's smoothing term
unchanged while changing the other candidate classes' terms. Thirteen focused
tests pass, including source guard and invalid-input checks.

## Scope limits

This validates only one candidate-class smoothing contribution. It does not
port the complete update loop, cached membership update order, inverse-energy
normalization, atlas shape priors, topology gate, deformation, or final
segmentation. The Java output is from a source-formula excerpt and is not full
algorithm equivalence evidence.

## Full-suite evidence

- Command: `./.venv/bin/pytest -q`
- Result: **291 passed in 53.02 seconds**
- Pytest: 9.1.1; Python: 3.13.2
- Tested implementation commit: `d8a02be707cc0a14ed1876e6bb7b3476aba71a43`
- Run date: 2026-09-28
