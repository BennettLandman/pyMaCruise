# Folded phantom: bounded continuing TOADS slice

This fixture is a deterministic, asymmetric folded-band phantom. The checked-in
NIfTI-1 `.nii.gz` images contain the starting segmentation, one-channel input,
initial class memberships, and transformed class shape values. Their array
axes are x/y/z in RAS order and the affine uses positive R/A/S axis directions
with spacings `(1.0, 1.25, 1.5)` mm. `generate_fixture.py` recreates the four
inputs.

The test composes the selected source order for one continuing outer-loop
iteration:

1. `MultipleToadSegmentation.propagateThinning()`;
2. one supported ordinary-OBJECT membership sweep using the carried closest
   label field from before this iteration;
3. `AlgorithmMultipleToads` float32 energy-ratio and non-stop decision;
4. `MultipleToadSegmentation.propagateGrowing()`;
5. relation, boundary-distance, and closest-label refresh from the grown labels.

The full-output field snapshot is in `composition_snapshot.json`. It is an
exact Python regression snapshot of the integrated call (including float32
bit patterns and every label/membership/relation voxel); it is not an integrated
Java output or a claim of full TOADS equivalence. The component algorithms
already have separate selected-Java comparisons in
`test_toads_initialized_thinning_pass.py`,
`test_toads_object_membership_pass.py`,
`test_toads_iteration_control.py`, `test_toads_growing_pass.py`, and the
relation tests.

The selected upstream source files are from VandyCRUISE:

- `VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmMultipleToads.java`
  SHA-256 `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`
- `VandyCRUISE/src/edu/jhmi/rad/medic/methods/MultipleToadSegmentation.java`
  SHA-256 `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`

This is deliberately only an integrated iteration slice. It accepts the
ordinary OBJECT model parameters rather than porting the upstream centroid and
variance dispatcher, excludes atlas registration and optional inhomogeneity
correction, and does not establish whole-run stopping/convergence or clinical
performance. Lesion/outlier membership dispatch is also outside this slice.
