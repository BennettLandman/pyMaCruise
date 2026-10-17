# AtlasEM descriptor-loaded rigid registration

## Scope

Commit `9a3a61a79c1523093c9e28999f5c22ef1b763cc9` adds `register_atlas_em_descriptor_pyramid`. It parses the selected descriptor, loads topology and class shape maps from `.nii.gz`, assembles source class order, passes the registered-class flags to the selected rigid registration adapter, and returns the descriptor-loaded atlas with the registration history. Subject memberships, image spacing, pyramid level count, initial transform, and iteration cap remain explicit inputs; the template spacing defaults to the descriptor's shape spacing.

This composes the Python atlas-loading and registration slices into an executable synthetic stage. It does not initialize subject memberships, run EM, execute the following Dura step, or establish Java optimizer parity. The missing `ibsr8obj` atlas payload remains a prerequisite for the selected real configuration.

## Validation

The focused descriptor-to-registration test passed **1 test in 0.58 s** and confirms that one scale-2 registration iteration ran with the two classes and registered-class vector in descriptor order. The broader AtlasEM regression passed **34 tests, 3 expected orientation warnings, in 6.06 s** on tested commit `9a3a61a79c1523093c9e28999f5c22ef1b763cc9`. This uses committed test code but only synthetic data. No full project suite or direct Java output comparison was run.
