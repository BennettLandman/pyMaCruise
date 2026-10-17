# AtlasEM descriptor, registration, and native recrop handoff

## Scope

Commit `9abb80e88d720d3261e9a643c5e7faadcad61b25` adds `register_and_recrop_atlas_em_from_descriptor`. It loads the `.nii.gz` atlas from the selected text descriptor, preserves descriptor class order and registration flags, runs the configured rigid pyramid adapter, then passes its final transform through the existing source-shaped uncrop, foreground-bound recomputation, and recrop sequence.

The adapter consumes subject memberships and normalized native-grid images prepared by earlier stages. It does not implement AtlasEM subject initialization or EM updates and has not been compared to a direct execution of `DemonToadDeformableAtlas`. The synthetic integration completes one scheduled optimizer iteration and the native-grid recrop only. Selected atlas assets remain unavailable.

## Validation

The focused descriptor-registration tests passed **2 tests in 0.60 s**. The broader AtlasEM regression across descriptor loading, registration and recropping, statistics, memberships, preprocessing, outputs, CLI handoff, and synthetic Dura integration passed **35 tests, 3 expected orientation warnings, in 6.03 s** on tested commit `9abb80e88d720d3261e9a643c5e7faadcad61b25`. No full suite or direct Java optimizer comparison was run.
