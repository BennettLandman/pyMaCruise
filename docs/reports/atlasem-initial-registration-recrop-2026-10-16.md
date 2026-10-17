# AtlasEM initial membership, registration, and recrop integration

## Scope

Commit `2010e4cecaf868fd95815f1b4027dd201d6506a2` adds `initialize_register_recrop_atlas_em_from_descriptor`. Given a prepared image crop, its native-grid normalized image data and robust bounds, a descriptor, and source-order modality names, it computes descriptor-driven initial memberships, loads class-ordered RAS shape priors, runs the selected rigid pyramid adapter, and recrops the image on its native grid using the registered transform.

This creates a restartable stage boundary from preprocessed/cropped input through initial memberships, one or more configured rigid updates, and post-registration recropping. It does not prepare or read the subject MRI, derive the user brain mask, run the full AtlasEM EM loop, generate the full output set, or feed the Dura stage. The synthetic descriptor does not represent the unavailable selected `ibsr8obj` atlas. No direct Java registration equivalence is claimed.

## Validation

The focused descriptor initialization/registration tests passed **3 tests in 0.63 s**. The broader AtlasEM regression passed **36 tests, 3 expected orientation warnings, in 6.01 s** on tested commit `2010e4cecaf868fd95815f1b4027dd201d6506a2`. The integration test asserts descriptor class ordering, initialized membership shape and finiteness, one scheduled scale-2 registration update, and finite native-grid recrop output. No full project suite was run.
