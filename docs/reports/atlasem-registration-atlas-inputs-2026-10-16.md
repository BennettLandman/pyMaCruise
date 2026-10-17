# AtlasEM descriptor-to-registration input assembly

## Scope

Commit `a8b48b6e5910a66ed33d9db923a917ad61a2f200` adds `load_atlas_em_registration_inputs`. It loads the descriptor's `.nii.gz` topology and class shape maps through the RAS loader, orders shape priors by the descriptor's class list, and returns the stack in `(x, y, z, class)` order together with the source class names, labels, and registered-shape flags. It rejects incomplete or unexpected class shape maps rather than assigning maps to the wrong class.

The result matches the array organization expected by the existing rigid registration adapter. It does not initialize subject memberships, run the optimizer, or validate the source-selected `ibsr8obj` payload.

## Validation

The descriptor-focused synthetic suite passed **10 tests in 0.63 s**. It checks exact class order, labels, flags, and the channel-last stack values after a synthetic LAS-to-RAS axis flip. The full focused AtlasEM group passed **33 tests, 3 expected orientation warnings, in 6.04 s** on tested commit `a8b48b6e5910a66ed33d9db923a917ad61a2f200`. No full project suite or direct Java optimizer comparison was run.
