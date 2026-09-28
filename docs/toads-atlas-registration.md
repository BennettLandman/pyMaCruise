# TOADS atlas registration port notes

**Status:** Investigation started 2026-09-28. This document identifies the
selected rigid registration behavior and the next source-comparison work. It
does not claim that atlas registration is implemented.

## Selected execution path

The selected VandyCRUISE `AlgorithmMultipleToads` path computes exact tissue
memberships (`useApprox=false`), initializes the atlas transform from
`MultipleToadPreprocess.getTransform()`, calls
`MultipleToadDeformableAtlas.initShapeRegistration(..., "rigid", ...)`, runs
`updateShapeRegistration(...)`, then calls `registerShapesPyramid()`. The
driver subsequently recalculates transformed support bounds and repeats the
crop-boundary calculation. In the selected wrapper configuration, the
registration transform is rigid. This is the atlas-to-native-input alignment
inside TOADS; the user's decision to work in native T1 RAS space excludes the
historical SPECTRE affine preprocessing.

## Objective and update rule

For exact memberships, `computeExactRegistrationCoefficients` scans the
selected subsampled input grid and each registered non-background class. At a
mapped coordinate it samples a shape prior `p_k` with closest-edge linear
interpolation and skips the contribution when `p_k <= 0`. The source then
accumulates

```text
cost += membership_k^2 * p_k^2
norm += p_k^2
energy = cost / norm
```

The optimizer maximizes this normalized overlap-like energy. It builds a
diagonal gradient/Hessian approximation from the sampled prior's spatial
derivatives and rigid-transform coordinate derivatives. Each trial updates
parameters by `trial[n] = transform[n] + lambda / hessian[n] * gradient[n]`.
An energy-improving trial is accepted and multiplies `lambda` by 1.5; a rejected
trial divides `lambda` by 1.5. Search stops after more than ten rejected
trials. The outer loop is bounded by the configured iteration count and stops
when relative energy change is at most `1e-6` or `lambda <= 0.001`. The source
uses float32 arithmetic and source-specific rotation/interpolation routines.

### Implemented objective-kernel slice

`src/pymacruise/toads_registration.py::exact_membership_registration_energy`
now implements that energy accumulation for an already sampled prior field,
including the source's first-class skip, registered-class selection, positive
prior filter, offset/subsample traversal, float32 accumulation, and empty-norm
NaN behavior. Two deterministic Java harness cases cover all registered
classes and a selected-class/subsampled scan. The harness is a source-formula
excerpt: it does not execute the full private source method, transform mapping,
interpolation, derivatives, or the rigid update.

## Python-library reuse decision

ANTsPy registration is already available for workflow stages whose source
contract is image-to-image registration. It is not a direct substitute here:
the TOADS objective is membership-weighted atlas-prior overlap with a custom
parameter update, rather than an ANTs image similarity metric and optimizer.
Using ANTs for this internal alignment would change the algorithm and requires
end-to-end evidence before it could be accepted as a substitute.

SciPy `ndimage.map_coordinates(order=1, mode="nearest", prefilter=False)` was
compared against the existing source-ordered float32 interpolation kernel on
10,000 deterministic random coordinates over a seeded 7×6×5 float32 field.
The largest observed absolute difference was one float32 ULP
(`1.1920928955078125e-7`), affecting 2,810 samples. This is a close numerical
match, but not exact. Since the atlas prior contributes repeatedly to
registration and segmentation, the current exact-order kernel remains the
reference path until a broader synthetic pipeline comparison shows the SciPy
difference stays within the project acceptance criteria. This avoids
reimplementing functionality where a library already matches while also
preserving source fidelity where it does not.

## Next implementation slices

1. Make a small Java harness for the exact rigid objective's energy, gradient,
   and diagonal Hessian on analytic synthetic membership/prior fields. Record
   transform, interpolation, derivative, accumulation order, and source hashes.
2. Implement and compare the objective/coefficients in Python, using SciPy
   only for operations shown to preserve the reference within measured limits.
3. Compare one rigid update and its accept/reject schedule before adding the
   full iteration loop or pyramid levels.
4. Validate transformed-bound refresh and registration-driven crop changes on
   larger asymmetric folded synthetic phantoms. Resolve the atlas descriptor
   versus NIfTI spacing discrepancy before any atlas execution claim.

## Source provenance

Selected source: `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java`,
SHA-256 `e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`.
Driver: `VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java`.
The exact local source revisions and wrapper parameter mapping remain recorded
in [toads-port-plan.md](toads-port-plan.md) and the source asset manifest.
