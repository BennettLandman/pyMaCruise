# Full-suite evidence — P18ck TOADS registration pyramid

- Date: 2026-09-30
- Tested implementation commit: `39035a15ec342624090df6bd7ebf15f76460b5e2`
- Command: `.venv/bin/pytest -q`
- Result: **629 passed**, 2 expected RAS reorientation warnings, in 129.93 seconds.
- Focused validation: 19 registration-sampling tests passed; the new selected-Java pyramid comparison passed exactly.
- Java reference regeneration: `sh scripts/generate_java_toads_registration_pyramid_reference.sh`.

The Java harness reads the same committed 17×15×13 membership and 19×17×15
atlas-prior `.nii.gz` fixtures as Python. The grids have three classes, odd
dimensions, and anisotropic RAS spacing. The reference uses the preserved
VandyCRUISE `ImageFunctions.subsample`, `ParametricTransform`,
`RotationMatrix`, interpolation, and derivative methods, with the selected
TOADS coefficient normalization and line-search control flow transcribed in
the harness.

At scales 4 and 2, Python exactly matches the harness's coarse dimensions,
starting transforms, final transforms, float32 energy/lambda/change values,
accepted flags, rejection counts, iteration counts, and stop decisions. Scale
4 accepts two updates and hands the resulting transform to scale 2; scale 2's
selected line search rejects all 11 trials, matching the source rejection
condition. This is evidence for the coarse-to-fine rigid-registration slice
only; it does not establish atlas physical alignment, full TOADS segmentation,
or full MaCRUISE equivalence.

Both warnings are from established tests that intentionally supply LPS
topology volume and paint-mask inputs, which are reoriented to RAS.
