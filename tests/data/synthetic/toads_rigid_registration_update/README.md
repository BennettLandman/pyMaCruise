# TOADS rigid registration update reference

The 4×4×4 membership grid and 5×5×5 atlas fields use deterministic linear and
cross-term formulas. Two registered foreground classes contribute to the
source normalized-overlap energy. The test exercises one improving line-search
update and one update that rejects all eleven trial steps and retains the
initial transform. A third case exercises three consecutive selected
`registerAllShapes` outer-loop updates: two accepted transforms followed by an
11-trial rejection. Rotations and translations are nonzero for the backtracked
case.

The test-only Java harness calls the supplied `RotationMatrix`,
`ParametricTransform`, and `ImageFunctions` implementations for coordinate
mapping, derivatives, and interpolation, and reproduces the selected exact
membership coefficient loop and `registerGradientDescent` call order. It
records initial energies/coefficient arrays, each trial transform and energy,
and final optimizer state after each outer-loop iteration. Source hashes for
the selected atlas class and the three copied Java helper classes are embedded
in `java_reference.json`.
Source-faithful coefficient normalization retains
the selected source's observed gradient/hessian call-order behavior. No atlas
or subject image data are included.

Regenerate with
`scripts/generate_java_toads_rigid_registration_update_reference.sh`. The
outer-loop case is a selected source-formula harness, not execution of the
whole VandyCRUISE application. These checks do not establish full atlas
alignment, pyramid optimization, deformable TOADS, or segmentation
equivalence.
