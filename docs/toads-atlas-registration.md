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

### Coefficient normalization call-order finding

The final exact-registration coefficient loop in the selected source first
replaces each `hessian[i]` with `registrationMeasureHessian(hessian[i],
nhess[i], norm, normalize)`, whose active implementation returns `hess / norm`.
It then assigns `gradient[i] = registrationMeasureGradient(hessian[i],
nhess[i], norm, normalize)` rather than passing the accumulated `gradient[i]`.
The active gradient helper returns its first argument divided by `norm`, so the
observed result is the already-scaled Hessian divided by `norm` again. The
source-faithful option in `exact_registration_measure_coefficients` preserves
that call order by default. An apparent intended `gradient / norm` result is
available only with `source_faithful=False`; no evidence yet supports using
that corrected calculation for historical output reproduction.

Two fixed coefficient vectors match a separate Java method-body excerpt
exactly after float32 conversion. This establishes the selected method's
call-order behavior, not the coefficient accumulators, optimizer schedule, or
alignment quality.

### Trilinear shape-prior derivatives

`trilinear_prior_derivatives` now follows the three selected
`ImageFunctions.linearInterpolation*derivative` methods. Its edge behavior is
not the energy sampler's closest-edge behavior: if any coordinate lies below
zero or above `size - 2`, all three derivatives are zero. Four deterministic
points match the source-formula Java excerpt exactly after float32 conversion,
including the upper face and both outside-bound cases. The rigid transform
Jacobian is now implemented as the next separate primitive.

### Rigid transform coordinate Jacobian

`rigid_image_to_template_jacobian` combines the selected rotation-matrix
derivatives with the RIGID branch of `ParametricTransform` for one voxel. It
returns a 3×6 float32 matrix (template x/y/z by three rotation and three
translation parameters). A fixed non-singular synthetic rotation and
anisotropic image/template spacing match the Java formula excerpt exactly,
including translation derivatives. The singular unit-norm rotation case is
rejected explicitly; it has not been shown to occur in the selected optimizer.
This does not yet validate transform mapping, prior-gradient composition, or
registration coefficient accumulation.

### Per-voxel exact coefficient accumulation

`exact_registration_coefficient_contributions` composes each sampled prior's
x/y/z derivative with the 3×6 rigid coordinate Jacobian, then accumulates the
selected source's cost, norm, six-entry gradient and diagonal Hessian, plus
norm-gradient/norm-Hessian terms. A deterministic four-class voxel case
matches all returned float32 accumulators from the Java method-formula
excerpt exactly. The background, an unregistered class, and a zero-prior class
exercise source skip rules. This validates one voxel's coefficient arithmetic
only; it does not validate full-volume loop ordering, transform-dependent
prior sampling, source coefficient normalization, line search, or alignment.

### One source registration update and backtracking step

`exact_registration_gradient_descent_step` applies the source diagonal update
and maximization rule: strict improvement accepts and grows lambda by 1.5;
otherwise lambda shrinks by 1.5, with termination after the 11th consecutive
rejection under the selected `itSupp=10` setting. Synthetic accepted and
rejected energy callbacks match the Java formula excerpt for parameter state,
energy, lambda, relative change, and reject count. This implements one
registration update around a supplied energy callback; it does not yet run the
TOADS energy or coefficient calculation end to end.

### Full sampled-grid coefficient accumulation

`exact_membership_registration_coefficients` now runs the sampled-volume
x/y/z/class accumulation in the selected Java loop order. It accepts arrays of
already sampled priors, their derivatives, and per-voxel rigid Jacobians, then
returns source-order float32 cost, norm, and coefficient accumulators. Two
small synthetic grids (full scan and offset/subsample scan) match the
standalone Java source-formula harness exactly across all returned values.
This is a numerical loop kernel only: it does not build transform-dependent
prior fields or execute the registration optimization. The Python nested loops
are currently correctness-first; runtime performance on full clinical grids
has not been measured.

### Integrated optimizer and streamed evaluation

`register_toads_atlas_rigid` connects the transformed class-prior sampler,
exact-membership energy and coefficient accumulation, source-faithful
normalization, and rigid line search. `source_faithful=True` remains the
default. The deterministic synthetic test repeats the loop and compares all
history fields exactly, exercises accepted trials and the iteration, energy,
and lambda stopping gates, and rejects an objective with no positive prior
support. This is optimizer wiring evidence only: no full Java execution of
`MultipleToadDeformableAtlas` or known-transform alignment claim is available.

P18u initially used dense arrays. P18v switched the optimizer to
`stream_toads_registration_energy` and
`stream_toads_registration_coefficients`, which process one voxel and class
at a time and keep source x/y/z/class accumulation order. The streamed path
matches dense energy, cost, normalization, gradient, Hessian, norm-gradient,
and norm-Hessian exactly on the anisotropic synthetic case, for both full and
offset/subsample scans. The optimizer's deterministic/stopping tests run over
the streaming implementation.

The synthetic benchmark measures one coefficient plus energy scan for three
classes: 0.1053 s at 6³, 0.2481 s at 8³, and 0.4903 s at 10³. The separate
Python-traced temporary allocation peak stays near 0.006 MiB and excludes the
preallocated inputs. Those are P18v scalar-path measurements, retained for
comparison; the current optimizer uses the P18w vectorized slab implementation.
It matches the scalar and dense paths exactly on two synthetic anisotropic
cases and takes 0.00689 s at 6³, 0.01176 s at 10³, and 0.02076 s at 16³ for
one coefficient plus energy scan (15.05×, 41.02×, and 97.83× speedups). Its
traced scratch peaks range from 0.0323 to 0.1580 MiB over those grids. These
small-grid results still do not establish full atlas runtime or full-run
performance. Details and the command are in [the P18w benchmark report](reports/toads-registration-vectorization-2026-09-28.md).

P18x adds the selected `registerShapesExactPyramid` schedule and reproduces
`ImageFunctions.subsample` block means, including float32 offset-loop order and
discarded incomplete blocks; its 9×8×7×3 asymmetric fixture matches the Java
formula harness exactly at factors 2 and 4. P18y extends rigid prior sampling
to the source `ParametricTransform.imageToTemplate` scale formula and checks
the existing scaled Jacobian against a Java formula harness at scales 1, 2, and
4. These checks use odd full-resolution dimensions, explicit full-grid image
and template centers, anisotropic spacing, and a nonzero rigid transform. They
validate the coordinate and Jacobian primitives only; multi-resolution
optimizer wiring and full Java optimizer equivalence remain open.

P18z wires the source levels into
`register_toads_atlas_rigid_pyramid`: each input is block-mean downsampled
from its native grid, the full-resolution centers are retained, scaled maps
and Jacobians drive the existing optimizer, and each level's final transform
initializes the next level. It resets the line-search step at each level,
matching `registerShapesExactPyramid`. A deterministic 17×15×13 membership /
19×17×15 atlas fixture runs scales 4 and 2, accepts updates at the first
level, checks exact transform handoff, and reproduces the full result exactly
on repeat. `levels=1` returns the input transform because the selected source
does not optimize a scale-1 stage. A separate scale-2 kernel comparison against
scalar mapped priors, derivatives, and Jacobians has exact energy agreement;
raw coefficient arrays are within `rtol=2e-7`, `atol=3e-10`, with maximum
absolute difference `1.49e-8` (one float32 ULP at the affected magnitude).
This is synthetic optimizer-integration evidence only, not a Java full-method
registration comparison or alignment-quality claim.

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

### Remaining implementation and validation slices

The current Python loop is integrated only on small synthetic arrays. Continue
with these source-grounded items before using bundled atlas grids or proceeding
to full segmentation:

1. Build a Java harness that executes the selected full rigid-registration
   method body and compare transformed coordinates, energy, coefficients,
   accepted/rejected trials, and the final transform on known asymmetric
   synthetic alignment cases. The existing harnesses compare isolated formula
   slices and the Python multi-level integration, not the private Java atlas
   registration method.
2. Exercise transformed support-bound refresh and registration-driven crop
   changes on larger asymmetric folded synthetic phantoms.
3. Resolve the descriptor (0.9375 mm) versus bundled image metadata (1.0 mm)
   geometry contract in a Java/Python atlas loader comparison before any atlas
   execution claim.
4. Integrate rigid registration into TOADS preprocessing, then continue the
   deformable membership/topology evolution and inter-pass dura-removal port;
   this document does not close P18 or unblock P34.

The rigid optimizer currently preserves the source coefficient normalization
call order by default. It does not establish that the source call order is
numerically stable; the opt-in corrected calculation remains experimental and
must stay out of defaults without owner sign-off.

## Source provenance

Selected source: `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java`,
SHA-256 `e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2`.
Driver: `VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java`.
The exact local source revisions and wrapper parameter mapping remain recorded
in [toads-port-plan.md](toads-port-plan.md) and the source asset manifest.
