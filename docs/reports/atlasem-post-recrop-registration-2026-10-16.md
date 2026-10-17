# AtlasEM post-recrop single-scale registration

## Scope

This slice implements the second rigid-registration call in the selected
`AlgorithmAtlasEMSegmentation` flow. After AtlasEM initializes priors and
recomputes subject memberships on the recropped grid, the Java wrapper calls
`updateShapeRegistration(..., mainRegIter, 1, ...)` followed by
`atlas.registerShapes()`.

## Source evidence

The selected VandyCRUISE source is
`toads-cruise/src/edu/jhmi/rad/medic/algorithms/AlgorithmAtlasEMSegmentation.java`
(the `initialize()` call sequence near lines 397–405) and
`toads-cruise/src/edu/jhmi/rad/medic/methods/DemonToadDeformableAtlas.java`.
The latter sets voxel stride to 2 in `updateShapeRegistration`, loops up to
`itMax` in `registerAllShapes`, and performs one exact coefficient calculation
per gradient-descent iteration. In the exact single-scale coefficient method,
the source first normalizes `hessian[i]`, then passes that same array element
as the gradient argument to `registrationMeasureGradient`. This historical
call order is preserved with `source_faithful=True`; the pyramid call has a
different source formula and remains a separate adapter.

The MATLAB wrapper supplies the maximum EM iteration count and registration
mode but does not override the Java plugin's registration settings. The
selected plugin defaults are 4 pyramid levels, 50 initial-registration
iterations, and 1 post-recrop iteration; see the [selected-defaults report](atlasem-selected-registration-defaults-2026-10-16.md).
The lower-level adapter still accepts iteration count, cropped memberships,
class bounds, geometry and current transform explicitly.

## Implementation and validation

`register_atlas_em_rigid_single_scale` in
`src/pymacruise/atlas_em_registration.py` connects the already source-shaped
rigid-registration kernel to the selected one-level AtlasEM stage.
`initialize_register_recrop_register_atlas_em_from_descriptor` composes the
initial descriptor memberships, pyramid update, native-grid recrop, second
membership initialization, transformed class bounds, and post-recrop update.
The synthetic descriptor test now exercises both registration stages on a
12×12×12 native grid, with one optimizer iteration at each scheduled stage,
finite memberships and final transform. A separate 7×7×7 optimizer case
exercises source-selected stride two and class registration filtering. The source
coefficient-call-order behavior is implemented and independently covered by
the existing selected Java rigid-update fixture. This report does not claim a
direct Java comparison of the `DemonToadDeformableAtlas` optimizer output.

Focused command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/pytest -q tests/test_atlas_em_registration_from_descriptor.py tests/test_toads_rigid_registration_update.py
```

Result: **5 passed in 0.62 s**, integration implementation commit
`db1a68c9364196115280df013b8b9ddba7558dfe` (the standalone adapter was first
introduced in `5a632d12fdf2e1791d34b477adbcbf930c63b05f`).

The synthetic case rejects its trial steps under the source-faithful
coefficient behavior and returns the unchanged finite transform. This is
recorded as an observed source-faithful synthetic outcome, not as evidence
that the selected atlas optimizer has converged. A direct Java harness for the
`DemonToadDeformableAtlas` post-recrop output comparison remains open, as do
the missing `ibsr8obj` files and full AtlasEM/Dura workflow integration.
